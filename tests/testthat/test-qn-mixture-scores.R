# The fit_mixture() finish takes each case's score from one E-step
# (R/qn_scores.R) instead of differencing the likelihood 2p times. Those scores
# are only worth having if they are the derivatives of the objective the
# finish climbs, in the order it packs its coordinates, so every family is
# held here to central differences of that same objective, to 1e-5, at the
# fit and at a point off it. A family without scores must say so (NULL) rather
# than hand the finish a wrong or misaligned matrix.

.sc_q <- function(e) suppressMessages(suppressWarnings(e))

.sc_bin <- function(n = 300, J = 4, K = 3, na = FALSE, seed = 1) {
  set.seed(seed)
  cl <- sample(seq_len(K), n, TRUE)
  P  <- matrix(seq(.15, .85, length.out = K), K, J)
  P[, 2] <- rev(P[, 2])
  X  <- t(vapply(cl, function(k) rbinom(J, 1, P[k, ]), numeric(J)))
  if (na) X[sample(length(X), 30)] <- NA
  colnames(X) <- paste0("u", seq_len(J))
  X
}
.sc_cont <- function(n = 300, J = 3, na = FALSE, seed = 2) {
  set.seed(seed)
  cl <- sample(1:3, n, TRUE)
  X  <- matrix(rnorm(n * J, sd = c(1, 1.5, .7)[cl]), n) + 2 * cl
  if (na) X[sample(length(X), 30)] <- NA
  X
}
.sc_count <- function(n = 300, J = 3, na = FALSE, seed = 3) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE)
  X  <- matrix(rpois(n * J, ifelse(cl == 1, 0.5, 3)), n)
  if (na) X[sample(length(X), 30)] <- NA
  X
}
.sc_cat <- function(n = 300, J = 3, na = FALSE, seed = 4) {
  set.seed(seed)
  cl <- sample(1:3, n, TRUE)
  pr <- list(c(.6, .2, .1, .1), c(.1, .6, .2, .1), c(.1, .1, .2, .6))
  X  <- t(vapply(cl, function(k) vapply(seq_len(J), function(j)
    sample.int(4, 1, prob = pr[[k]]), integer(1)), integer(J)))
  X[, J][X[, J] == 4] <- 3                  # a ragged item: three categories
  if (na) X[sample(length(X), 30)] <- NA
  X
}
.sc_fm <- function(X, K = 3, ...) .sc_q(fit_mixture(X, n_classes = K, n_init = 2,
                                                    random_state = 1, n_cores = 1, ...))
.sc_grp <- function(n = 300) factor(rep(1:3, length.out = n))

# Analytic against central differences, over the coordinates the finish moves.
# `prior = FALSE`: a family whose prior gradient is not written out, so the
# finish differences the prior; its builder must then say so (NULL, now or
# when evaluated, which is what .qn_finish_mixture() checks).
.sc_check <- function(fit, label, prior = TRUE) {
  Y  <- fit$Y
  pr <- .qn_mixture_problem(fit, fit$data, Y)
  expect_false(is.null(pr), label = label)
  free <- setdiff(seq_along(pr$par), pr$fixed)
  set.seed(11)
  pts <- list(fit = pr$par,
              off = pr$par + rnorm(length(pr$par), sd = 0.05))
  for (at in names(pts)) {
    par <- pts[[at]]
    an  <- .qn_mixture_scores(pr$unpack(par), pr$lay, fit$data, pr$Yp)
    lab <- paste(label, at)
    expect_false(is.null(an), label = lab)
    if (is.null(an)) next
    expect_identical(an$ll, pr$case_ll(par), label = paste(lab, "case LL"))
    fd <- .qn_fd_scores(pr$case_ll, par, free)
    expect_lt(max(abs(an$S[, free, drop = FALSE] - fd)), 1e-5, label = lab)
    # The prior's gradient and diagonal curvature. Second differences of a
    # prior of size ~100 are good to ~1e-4, so the curvature is held to 1e-3.
    pg <- .qn_mixture_prior_grad(fit, pr$lay, fit$data, Y, pr$wt)
    if (!prior) {
      expect_true(is.null(pg) || is.null(pg(pr$unpack(par))), label = lab)
      next
    }
    # A model whose prior is not priced (mixed items across groups) is never
    # finished, and has no prior gradient either.
    if (is.na(pr$prior(par))) { expect_null(pg, label = lab); next }
    expect_false(is.null(pg), label = paste(lab, "prior"))
    if (is.null(pg)) next
    a  <- pg(pr$unpack(par))
    nd <- .qn_fd_grad2(pr$prior, par)
    expect_lt(max(abs(a$g[free] - nd$g[free])), 1e-5, label = paste(lab, "prior gradient"))
    expect_lt(max(abs(a$h[free] - nd$h[free]) / pmax(1, abs(nd$h[free]))), 1e-3,
              label = paste(lab, "prior curvature"))
  }
}

test_that("measurement-model scores equal central differences", {
  .sc_check(.sc_fm(.sc_bin(), measurement = "binary"), label = "binary")
  .sc_check(.sc_fm(.sc_bin(na = TRUE), measurement = "binary"), label = "binary NA")
  .sc_check(.sc_fm(.sc_cat(), measurement = "categorical"), label = "nominal")
  .sc_check(.sc_fm(.sc_cat(na = TRUE), measurement = "categorical"),
            label = "nominal NA")
  .sc_check(.sc_fm(.sc_cont(), measurement = "gaussian"), label = "gaussian unit")
  .sc_check(.sc_fm(.sc_cont(na = TRUE), measurement = "gaussian"),
            label = "gaussian unit NA")
  .sc_check(.sc_fm(.sc_cont(), measurement = "continuous", variances_equal = FALSE),
            label = "gaussian free variances")
  .sc_check(.sc_fm(.sc_cont(), measurement = "continuous", variances_equal = TRUE),
            label = "gaussian equal variances")
  .sc_check(.sc_fm(.sc_cont(na = TRUE), measurement = "continuous",
                   variances_equal = FALSE), label = "gaussian NA")
  .sc_check(.sc_fm(.sc_count(), K = 2, measurement = "count"), label = "poisson")
  .sc_check(.sc_fm(.sc_count(na = TRUE), K = 2, measurement = "count"),
            label = "poisson NA")
  .sc_check(.sc_fm(cbind(.sc_bin(J = 3), .sc_cont(J = 2), .sc_cat(J = 2)),
                   measurement = list(binary = 1:3, continuous = 4:5,
                                      categorical = 6:7),
                   variances_equal = FALSE), label = "mixed")
  X <- .sc_bin(); d <- data.frame(z = rbinom(nrow(X), 1, .5), w = rnorm(nrow(X)) > 0)
  .sc_check(.sc_fm(X, measurement = "binary", predictors = ~ z, data = d,
                   predictors_items = list(u1 = ~ z, u3 = ~ z + w),
                   predictors_items_by_class = "u3"),
            label = "binary with item covariates")
  # Ordinal items with item covariates: a uniform slope, a slope by class and
  # a ragged three-category item, with and without missing cells.
  for (na in c(FALSE, TRUE)) {
    X <- .sc_cat(na = na); d <- data.frame(z = rbinom(nrow(X), 1, .5))
    colnames(X) <- paste0("o", 1:3)
    .sc_check(.sc_fm(X, measurement = "ordinal", cats = c(4L, 4L, 3L),
                     predictors = ~ z, data = d,
                     predictors_items = list(o1 = ~ z, o3 = ~ z),
                     predictors_items_by_class = "o3"),
              label = paste("ordinal with item covariates", if (na) "NA"),
              prior = FALSE)
  }
})

test_that("multiple-group measurement scores equal central differences", {
  g <- .sc_grp()
  .sc_check(.sc_fm(.sc_bin(), measurement = "binary", group = g,
                   group_effects = "both", n_steps = 1),
            label = "binary, groups both")
  .sc_check(.sc_fm(.sc_bin(na = TRUE), measurement = "binary", group = g,
                   group_effects = "measurement", group_invariant_items = c(1, 3)),
            label = "binary, invariant items")
  .sc_check(.sc_fm(cbind(.sc_bin(J = 3), .sc_cont(J = 2)), group = g,
                   measurement = list(binary = 1:3, continuous = 4:5),
                   group_effects = "measurement", group_invariant_items = c(2, 4),
                   variances_equal = FALSE), label = "mixed, invariant items")
  for (ip in list("means", "covariances", c("means", "covariances")))
    for (ve in c(TRUE, FALSE))
      .sc_check(.sc_fm(.sc_cont(), measurement = "continuous", group = g,
                       group_effects = "measurement", group_invariant_params = ip,
                       variances_equal = ve),
                label = paste("gaussian invariant", paste(ip, collapse = "+"),
                              "equal variances", ve))
})

test_that("structural-model scores equal central differences", {
  X <- .sc_bin(); n <- nrow(X); g <- .sc_grp(n)
  set.seed(21)
  z <- data.frame(z = rnorm(n)); zf <- data.frame(z = rnorm(n), f = rbinom(n, 1, .4))
  one <- function(fit, label) .sc_check(fit, label)
  one(.sc_fm(X, measurement = "binary", predictors = zf, n_steps = 1),
      "covariates")
  one(.sc_fm(X, measurement = "binary", group = g, group_effects = "prevalence",
             n_steps = 1), "prevalence free")
  one(.sc_fm(X, measurement = "binary", group = g, group_effects = "prevalence",
             group_prevalence_equal = 1, n_steps = 1), "prevalence, one class equal")
  one(.sc_fm(X, measurement = "binary", group = g, group_effects = "prevalence",
             group_prevalence_equal = 1:2, n_steps = 1), "prevalence, two classes equal")
  y <- rnorm(n) + (X[, 1] == 1); y[sample(n, 15)] <- NA
  one(.sc_fm(X, measurement = "binary", outcome = y, n_steps = 1), "distal continuous")
  one(.sc_fm(X, measurement = "binary", outcome = y, outcome_covariates = zf,
             n_steps = 1), "distal continuous, adjusted")
  one(.sc_fm(X, measurement = "binary", outcome = y, outcome_covariates = zf,
             slopes = "class_specific", n_steps = 1),
      "distal continuous, class-specific slopes")
  one(.sc_q(fit_mixture_internal(
        X, cbind(y = y, z = zf$z), n_components = 3, measurement = "binary",
        structural = "distal_continuous_regression", n_steps = 1,
        n_init = 2, random_state = 1, n_cores = 1)), "distal continuous regression")
  yc <- factor(sample(1:3, n, TRUE)); yc[sample(n, 15)] <- NA
  one(.sc_fm(X, measurement = "binary", outcome = yc, n_steps = 1), "distal categorical")
  one(.sc_fm(X, measurement = "binary", outcome = yc, outcome_covariates = z,
             n_steps = 1), "distal categorical, adjusted")
  one(.sc_fm(X, measurement = "binary", outcome = yc, outcome_covariates = z,
             slopes = "class_specific", n_steps = 1),
      "distal categorical, class-specific slopes")
})

test_that("a frozen measurement block is left out", {
  # The two-step estimator's step 2 (R/stepwise.R) holds the measurement block
  # where step 1 left it; its columns are zeros and are never read.
  X <- .sc_bin(); set.seed(22)
  f2 <- .sc_fm(X, measurement = "binary", predictors = rnorm(nrow(X)), n_steps = 1)
  f2$frozen <- "mm"
  .sc_check(f2, "measurement frozen")
  pr <- .qn_mixture_problem(f2, f2$data, f2$Y)
  expect_length(pr$fixed, length(pr$lay$idx$mm))
  S  <- .qn_mixture_scores(pr$unpack(pr$par), pr$lay, f2$data, pr$Yp)$S
  expect_true(all(S[, pr$fixed] == 0))
})

test_that("growth mixture scores equal central differences", {
  # Every Psi packing (shared, class variances with shared covariances, free,
  # a fixed element), residual variances by class and constant over
  # occasions, growth-factor regressions, a class regression and missing
  # occasions. The point off the fit moves the residual variances, which the
  # Psi unit is read from, so it also checks the unit stays the packing's.
  set.seed(6)
  n  <- 240
  Y  <- outer(rnorm(n, rep(c(0, 3), each = n / 2)), 0:3 * .5, "+") +
    rnorm(n, 0, .3) * outer(rep(1, n), 0:3) + matrix(rnorm(n * 4, sd = .5), n)
  Ym <- Y; Ym[sample(length(Ym), 60)] <- NA
  x  <- rnorm(n); z <- rbinom(n, 1, .5)
  fx <- data.frame(parameter = "growth_variance", item = "linear", class = 2,
                   value = 0.05)
  gm <- function(Y, ...) .sc_q(fit_gmm(Y, times = 4, n_classes = 2, n_init = 2,
                                       random_state = 1, n_cores = 1, ...))
  .sc_check(gm(Y), "growth, psi equal")
  .sc_check(gm(Ym, psi = "free"), "growth, psi free NA")
  .sc_check(gm(Y, psi = "variances", fix = fx), "growth, psi variances fixed")
  .sc_check(gm(Y, residual = "constant", residual_equal = FALSE),
            "growth, residual constant by class")
  .sc_check(gm(Y, growth_predictors = cbind(x, z), growth_predictors_equal = FALSE),
            "growth, regressions by class")
  .sc_check(gm(Y, predictors = data.frame(x = x), n_steps = 1),
            "growth, class regression")
})

test_that("growth mixture curvature equals differenced scores", {
  # The finish's closed-form Hessian (.qn_sn_hessian()) against central
  # differences of the analytic gradient, at the fit and off it, over the
  # packings it covers; then the finish's own hook, the class-size prior's
  # curvature included, and the models it leaves to the default schedule.
  set.seed(6)
  n  <- 240
  Y  <- outer(rnorm(n, rep(c(0, 3), each = n / 2)), 0:3 * .5, "+") +
    rnorm(n, 0, .3) * outer(rep(1, n), 0:3) + matrix(rnorm(n * 4, sd = .5), n)
  Ym <- Y; Ym[sample(length(Ym), 60)] <- NA
  fx <- data.frame(parameter = "growth_variance", item = "linear", class = 2,
                   value = 0.05)
  off <- list(latent = 0, categorical = 0, poisson = 0, variances = 0)
  gm <- function(Y, ...) .sc_q(fit_gmm(Y, times = 4, n_classes = 2, n_init = 2,
                                       random_state = 1, n_cores = 1, ...))
  grad <- function(pr, X, p) {
    st <- pr$unpack(p)
    e  <- .qn_mixture_estep(st, pr$lay, X, NULL)
    colSums(.qn_mixture_score_matrix(st, pr$lay, X, NULL, e$R))
  }
  fd_hess <- function(f, p) {
    H <- vapply(seq_along(p), function(i) {
      h <- 1e-5 * max(1, abs(p[i]))
      a <- p; a[i] <- a[i] + h; b <- p; b[i] <- b[i] - h
      (f(a) - f(b)) / (2 * h)
    }, numeric(length(p)))
    (H + t(H)) / 2
  }
  check <- function(fit, label) {
    X  <- fit$data
    pr <- .qn_mixture_problem(fit, X, NULL)
    set.seed(11)
    for (at in c("fit", "off")) {
      p  <- if (at == "fit") pr$par else pr$par + rnorm(length(pr$par), sd = 0.02)
      st <- pr$unpack(p)
      e  <- .qn_mixture_estep(st, pr$lay, X, NULL)
      hs <- .qn_sn_hessian(st$mm, X, e$R, rep(1, nrow(X)), st$weights)
      lab <- paste(label, at)
      expect_lt(max(abs(hs$S - .qn_mixture_score_matrix(st, pr$lay, X, NULL, e$R))),
                1e-10, label = paste(lab, "scores"))
      Hd <- fd_hess(function(q) grad(pr, X, q), p)
      expect_lt(max(abs(hs$H - Hd)) / max(abs(Hd)), 1e-7, label = lab)
    }
  }
  check(gm(Y, bayes_constants = off), "growth curvature, psi equal")
  check(gm(Ym, psi = "free", bayes_constants = off), "growth curvature, psi free NA")
  check(gm(Y, psi = "variances", fix = fx, bayes_constants = off),
        "growth curvature, psi variances fixed")
  check(gm(Ym, psi = "free", residual = "constant", residual_equal = FALSE,
           bayes_constants = off), "growth curvature, residual constant by class")

  # The hook, with the default class-size prior: the information it hands the
  # finish is minus the Hessian of the penalised objective.
  fit <- gm(Ym, psi = "free")
  X   <- fit$data
  pr  <- .qn_mixture_problem(fit, X, NULL)
  pg  <- .qn_mixture_prior_grad(fit, pr$lay, X, NULL, pr$wt)
  pgp <- function(p) pg(pr$unpack(p))
  expect_gt(.bayes_alpha(fit$mm, "latent"), 0)
  hook <- .qn_sn_hessian_hook(fit, pr$lay, X, pr, pgp)
  expect_false(is.null(hook))
  st <- pr$unpack(pr$par)
  e  <- .qn_mixture_estep(st, pr$lay, X, NULL)
  c0 <- c(list(st = st, R = e$R), hook$pass(list(st = st, R = e$R)))
  free <- seq_along(pr$par)
  I  <- hook$info(c0, free)
  Hd <- fd_hess(function(q) grad(pr, X, q) + pgp(q)$g, pr$par)
  expect_lt(max(abs(I + Hd)) / max(abs(Hd)), 1e-7)

  # Left to the default schedule: growth covariates, a class regression.
  x <- rnorm(n)
  expect_null(.qn_closed_curvature(gm(Y, psi = "free", growth_predictors = cbind(x)),
                                   NULL))
  fr <- gm(Y, psi = "free", predictors = data.frame(x = x), n_steps = 1)
  expect_null(.qn_closed_curvature(fr, fr$Y))
  expect_false(is.null(.qn_closed_curvature(fit, NULL)))
})
