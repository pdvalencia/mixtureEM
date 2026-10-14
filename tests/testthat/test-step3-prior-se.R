# The class regression is fitted with the latent prior, so its variance is the
# penalised objective's (the note above .step3_scores(), R/step3_variance.R).
# Checked here without any external number: the prior's derivatives against
# finite differences of the objective .fit_mnl() maximises, stationarity of the
# penalised score at every estimator that reports this variance, and the two
# identities that define the estimators.

.ps_sim <- function(n = 400, seed = 64) {
  set.seed(seed)
  z1  <- rnorm(n)
  z2  <- rbinom(n, 1, 0.4)
  eta <- cbind(0.3 + 0.8 * z1 - 0.5 * z2, -0.2 - 0.4 * z1 + 0.6 * z2, 0)
  cls <- apply(exp(eta) / rowSums(exp(eta)), 1,
               function(p) sample.int(3, 1, prob = p))
  pis <- rbind(rep(0.85, 6), c(0.85, 0.85, 0.85, 0.15, 0.15, 0.15), rep(0.15, 6))
  X   <- t(vapply(cls, function(k) rbinom(6, 1, pis[k, ]), numeric(6)))
  list(X = X, d = data.frame(z1 = z1, z2 = z2))
}

# The step-3 pieces of a covariate fit `f` built on step-1 fit `fit`, in the
# frame the variance code works in: classes in the step-1 order with the
# regression's anchor moved last. `idx` picks the free block of
# f$sm$parameters$V_robust, which is laid out in f's own class order.
.ps_frame <- function(fit, f, Z) {
  K  <- nrow(f$sm$parameters$beta)
  D  <- ncol(Z)
  p1 <- fit$mm$parameters$pis
  pf <- f$mm$parameters$pis
  m  <- vapply(seq_len(K), function(k)
    which.min(rowSums(abs(sweep(p1, 2, pf[k, ])))), integer(1))  # f class k = step-1 class m[k]
  B  <- f$sm$parameters$beta
  a  <- which(rowSums(abs(B)) == 0)
  ff <- c(setdiff(seq_len(K), a), a)                               # f classes, anchor last
  r  <- exp(fit$log_resp)[, m[ff]]
  Cn <- .step3_classification_table(r, rep(1, nrow(r)))
  list(B = B[ff, ], r = r, Cn = Cn,
       idx = as.vector(t(outer((ff[-K] - 1L) * D, seq_len(D), "+"))))
}

test_that("the prior's derivatives match finite differences of its log-density", {
  set.seed(1)
  Z    <- cbind(1, rnorm(30), rbinom(30, 1, 0.5))
  Z    <- rbind(Z, Z[1:5, ])                    # duplicated patterns count once
  beta <- rbind(c(0.2, -0.4, 0.7), c(-0.3, 0.5, 0.1), 0)
  a    <- 1.5
  K <- 3; U <- unique(Z); U0 <- nrow(U)
  as_b <- function(p) rbind(matrix(p, K - 1, 3, byrow = TRUE), 0)
  lp <- function(p) (a / (K * U0)) * sum(log(softmax_rows(U %*% t(as_b(p)))))
  p0 <- as.vector(t(beta[1:2, ]))
  d  <- .covariate_prior_derivs(beta, Z, a)
  h  <- 1e-5
  num_g <- vapply(seq_along(p0), function(m) {
    e <- replace(numeric(length(p0)), m, h); (lp(p0 + e) - lp(p0 - e)) / (2 * h)
  }, numeric(1))
  num_H <- vapply(seq_along(p0), function(m) {
    e <- replace(numeric(length(p0)), m, h)
    (.covariate_prior_derivs(as_b(p0 + e), Z, a)$gradient -
       .covariate_prior_derivs(as_b(p0 - e), Z, a)$gradient) / (2 * h)
  }, numeric(length(p0)))
  expect_equal(d$gradient, num_g, tolerance = 1e-7)
  expect_equal(d$hessian, num_H, tolerance = 1e-7)
  expect_true(all(.covariate_prior_derivs(beta, Z, 0)$hessian == 0))
})

test_that("a three-step fit ends at the penalised maximum and its variances are that objective's", {
  s   <- .ps_sim()
  fit <- fit_mixture(s$X, n_classes = 3, measurement = "binary", n_init = 5,
                     random_state = 1, n_cores = 1)
  Z   <- cbind(1, as.matrix(s$d))
  n   <- nrow(Z)
  w   <- rep(1, n)

  for (se in c("hessian", "robust")) {
    f  <- suppressMessages(add_covariates(fit, predictors = ~ z1 + z2,
                                          data = s$d, se = se))
    fr <- .ps_frame(fit, f, Z)
    pc <- .step3_pieces(fr$B, Z, fr$r, fr$Cn)
    pr <- .covariate_prior_derivs(fr$B, Z, 1)
    S  <- .step3_scores(pc, Z, w, 3)
    # Stationary on the penalised objective, and not on the plain one.
    expect_lt(max(abs(colSums(S) + pr$gradient)), 1e-6)
    expect_gt(max(abs(pr$gradient)), 1e-3)

    Hi <- solve(-(.step3_hessian(pc, Z, w, fr$r, 3) + pr$hessian))
    V  <- f$sm$parameters$V_robust[fr$idx, fr$idx]
    if (se == "hessian") {
      expect_equal(V, Hi, tolerance = 1e-6)
    } else {
      M <- n / (n - 1) * crossprod(sweep(S, 2, colMeans(S)))
      expect_equal(V, Hi %*% (M - pr$hessian) %*% Hi, tolerance = 1e-6)
    }
  }
})

test_that("without the latent prior the step-3 variance is the plain sandwich", {
  s   <- .ps_sim(seed = 65)
  fit <- fit_mixture(s$X, n_classes = 3, measurement = "binary", n_init = 5,
                     random_state = 1, n_cores = 1,
                     bayes_constants = list(latent = 0))
  f   <- suppressMessages(add_covariates(fit, predictors = ~ z1 + z2,
                                         data = s$d, se = "robust"))
  Z   <- cbind(1, as.matrix(s$d))
  n   <- nrow(Z)
  w   <- rep(1, n)
  fr  <- .ps_frame(fit, f, Z)
  pc  <- .step3_pieces(fr$B, Z, fr$r, fr$Cn)
  S   <- .step3_scores(pc, Z, w, 3)
  expect_lt(max(abs(colSums(S))), 1e-6)
  Hi  <- solve(-.step3_hessian(pc, Z, w, fr$r, 3))
  expect_equal(f$sm$parameters$V_robust[fr$idx, fr$idx],
               Hi %*% (n / (n - 1) * crossprod(sweep(S, 2, colMeans(S)))) %*% Hi,
               tolerance = 1e-6)
})

test_that("the two-step estimate is stationary on the penalised objective", {
  s   <- .ps_sim(seed = 66)
  f2  <- suppressMessages(fit_mixture(s$X, n_classes = 3, measurement = "binary",
                                      n_steps = 2, predictors = s$d, n_init = 5,
                                      random_state = 1, n_cores = 1))
  Y    <- as.matrix(s$d)
  info <- .twostep_information(f2, s$X, Y, rep(1, nrow(Y)))
  B    <- f2$sm$parameters$beta
  pg   <- .covariate_prior_derivs(sweep(B, 2, B[3, ], "-"),
                                  .covariate_design(f2$sm, Y), 1)$gradient
  expect_lt(max(abs(colSums(info$S2) + pg)), 1e-3)
  expect_gt(max(abs(colSums(info$S2))), 10 * max(abs(colSums(info$S2) + pg)))
})
