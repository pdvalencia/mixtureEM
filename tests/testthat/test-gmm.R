# ==============================================================================
# Growth mixture models (the structured multivariate-normal emission)
# ==============================================================================
#
# Checks that need no reference program: the likelihood against closed-form
# densities, boundary reporting, parameter-count arithmetic, and estimation
# on simulated data with known structure.

# ------------------------------------------------------------------------------
# Level 1: the likelihood
# ------------------------------------------------------------------------------

test_that("the structured-normal likelihood is the multivariate normal density", {
  set.seed(101)
  n <- 300; ts <- 0:3
  L <- .lcga_design(ts, 1)
  Y <- matrix(rnorm(n * 4, 2, 2), n, 4)

  alpha <- rbind(c(0.5, 0.30), c(2.0, -0.20))
  psi   <- matrix(c(0.8, 0.1, 0.1, 0.3), 2, 2)
  theta <- rbind(c(0.5, 0.6, 0.7, 0.8), c(0.5, 0.6, 0.7, 0.8))

  mm <- structured_normal_model(2, L, random_effects = "intercept_slope",
                                psi = "equal", residual = "occasion")
  mm$parameters$alpha <- alpha
  mm$parameters$psi   <- list(psi, psi)
  mm$parameters$theta <- theta

  # Written out the long way on purpose: a class contributes a single
  # multivariate normal density over the whole T-vector, whose covariance is
  # Lambda Psi Lambda' + Theta. That the occasions are *not* conditionally
  # independent given class is the entire difference from LCGA.
  expected <- matrix(0, n, 2)
  for (k in 1:2) {
    mu  <- as.vector(L %*% alpha[k, ])
    Sig <- L %*% psi %*% t(L) + diag(theta[k, ])
    ch  <- chol(Sig)
    z   <- backsolve(ch, t(sweep(Y, 2, mu, "-")), transpose = TRUE)
    expected[, k] <- -0.5 * (4 * log(2 * pi) + 2 * sum(log(diag(ch))) +
                               colSums(z^2))
  }

  expect_equal(log_likelihood(mm, Y), expected, tolerance = 1e-12)
})

test_that("with no random effects the emission is the Gaussian LCGA likelihood", {
  set.seed(103)
  n <- 200; ts <- 0:3
  L <- .lcga_design(ts, 1)
  Y <- matrix(rnorm(n * 4, 1, 1.5), n, 4)

  alpha <- rbind(c(0.5, 0.30), c(2.0, -0.20))
  disp  <- c(0.8, 1.5)

  sn <- structured_normal_model(2, L, random_effects = "none",
                                residual = "constant", residual_equal = FALSE)
  sn$parameters$alpha <- alpha
  sn$parameters$psi   <- rep(list(matrix(0, 0, 0)), 2)
  sn$parameters$theta <- cbind(disp, disp, disp, disp)

  lc <- lcga_model(2, design = L, family = "gaussian")
  lc$parameters$coefs      <- alpha
  lc$parameters$dispersion <- disp

  expect_equal(log_likelihood(sn, Y), log_likelihood(lc, Y), tolerance = 1e-10)
})

test_that("an unobserved occasion drops out of the likelihood", {
  set.seed(107)
  n <- 150
  L <- .lcga_design(0:3, 1)
  Y <- matrix(rnorm(n * 4), n, 4)
  Y_miss <- Y; Y_miss[1:20, 3] <- NA

  planted <- function(design) {
    mm <- structured_normal_model(2, design, random_effects = "intercept_slope",
                                  psi = "equal", residual = "occasion")
    mm$parameters$alpha <- rbind(c(0.3, 0.2), c(0.9, -0.1))
    psi <- matrix(c(0.7, 0.05, 0.05, 0.2), 2, 2)
    mm$parameters$psi   <- list(psi, psi)
    mm$parameters$theta <- rbind(c(0.4, 0.5, 0.6, 0.7), c(0.4, 0.5, 0.6, 0.7))
    mm
  }

  full <- log_likelihood(planted(L), Y_miss)

  # The same cases with occasion 3 removed from the model entirely must give the
  # same contribution. For a model with random effects this is a stronger claim
  # than it is for LCGA: the marginal of a multivariate normal over the observed
  # occasions has to equal the model built on those occasions alone, which is
  # true only because Sigma is being sub-setted rather than the density masked.
  mm3 <- planted(L[-3, , drop = FALSE])
  mm3$parameters$theta <- mm3$parameters$theta[, -3, drop = FALSE]
  drop <- log_likelihood(mm3, Y_miss[1:20, -3, drop = FALSE])

  expect_equal(full[1:20, ], drop, tolerance = 1e-12)
})

test_that("a solution on the boundary is reported rather than passed off", {
  L  <- .lcga_design(0:3, 1)
  mm <- structured_normal_model(2, L, psi = "free", residual = "occasion",
                                residual_equal = FALSE)
  mm$parameters$alpha <- rbind(c(1, 0.5), c(3, 0.1))
  mm$parameters$psi   <- list(matrix(c(1, 0.1, 0.1, 0.3), 2, 2),
                              matrix(c(1, 0.1, 0.1, 0.3), 2, 2))
  mm$parameters$theta <- rbind(c(0.5, 0.5, 0.5, 0.5), c(0.5, 0.5, 0.5, 0.5))
  expect_length(.gmm_boundary(mm), 0L)

  floored <- mm
  floored$parameters$theta[2, 3] <- .sn_theta_floor
  expect_match(.gmm_boundary(floored), "residual variance at zero")
  expect_match(.gmm_boundary(floored), "class 2, occasion 3")

  degenerate <- mm
  degenerate$parameters$psi[[1]] <- matrix(0, 2, 2)
  expect_match(.gmm_boundary(degenerate), "singular in class 1")
  expect_match(.gmm_boundary(degenerate), "fit_lcga")
})

test_that("fitting random effects to data that have none warns and collapses", {
  # No within-class variation in the growth factors: every case sits on its
  # class's line up to occasion noise, which is an LCGA. Psi has nowhere to go
  # but zero, and a user who asked for a growth mixture model should be told
  # that is what happened rather than reading a variance of 0.000 as an
  # estimate.
  set.seed(151)
  n <- 400; Tn <- 4
  cls <- rbinom(n, 1, 0.5)
  Y <- outer(ifelse(cls == 1, 3, 0), rep(1, Tn)) +
    outer(ifelse(cls == 1, 0.8, -0.2), 0:3) +
    matrix(rnorm(n * Tn), n, Tn)

  expect_warning(
    fit <- fit_gmm(Y, times = 4, n_classes = 2, n_init = 5, random_state = 3),
    "boundary of the parameter space")
  expect_match(fit$growth$boundary, "singular in class")
  expect_lt(max(diag(fit$growth$psi[[1]])), 0.05)
})

# ------------------------------------------------------------------------------
# Behaviour of the emission and the wrapper
# ------------------------------------------------------------------------------

# Two classes with genuine within-class variation in both growth factors.
#
# Simulating without random effects would be a trap rather than a simplification:
# Psi would converge to zero, which is a boundary solution the fit warns about,
# and every assertion that involves Psi would then be testing the degenerate
# case instead of the model.
.sim_gmm <- function(n = 600, Tn = 4, prob = 0.5, sd_i = 1.0, sd_s = 0.35,
                     sd_e = 0.7) {
  cls   <- rbinom(n, 1, prob)
  icept <- ifelse(cls == 1, 3, 0) + rnorm(n, 0, sd_i)
  slope <- ifelse(cls == 1, 0.8, -0.2) + rnorm(n, 0, sd_s)
  outer(icept, rep(1, Tn)) + outer(slope, seq.int(0L, Tn - 1L)) +
    matrix(rnorm(n * Tn, 0, sd_e), n, Tn)
}

test_that("free parameters are counted from the constraints, not the storage", {
  L <- .lcga_design(0:4, 1)
  np <- function(...) n_parameters(structured_normal_model(3, L, ...))

  # 3 classes x 2 growth means = 6, then the covariance structure.
  expect_equal(np(random_effects = "intercept_slope", psi = "equal",
                  residual = "constant", residual_equal = TRUE), 6 + 3 + 1)
  expect_equal(np(random_effects = "intercept_slope", psi = "equal",
                  residual = "occasion", residual_equal = TRUE), 6 + 3 + 5)
  expect_equal(np(random_effects = "intercept_slope", psi = "free",
                  residual = "occasion", residual_equal = FALSE),
               6 + 3 * 3 + 5 * 3)
  expect_equal(np(random_effects = "intercept", psi = "equal",
                  residual = "occasion", residual_equal = TRUE), 6 + 1 + 5)
  expect_equal(np(random_effects = "none", psi = "equal",
                  residual = "occasion", residual_equal = TRUE), 6 + 0 + 5)
})

test_that("the identification arithmetic is right, and fit_gmm() stays inside it", {
  # A class supplies T(T+1)/2 distinct covariance elements and spends
  # q(q+1)/2 + (1 or T) of them.
  expect_error(.gmm_check_identified(4L, 4L, FALSE), "not identified")
  expect_error(.gmm_check_identified(2L, 2L, FALSE), "not identified")
  expect_silent(.gmm_check_identified(4L, 3L, FALSE))   # 10 spent of 10
  expect_silent(.gmm_check_identified(4L, 2L, FALSE))   # 7 spent of 10
  expect_silent(.gmm_check_identified(4L, 2L, TRUE))    # 4 spent of 10

  # With a polynomial design the check can never actually fire, and it is worth
  # knowing why rather than discovering it as dead code: fit_gmm() already
  # requires degree + 1 < times, and every random effect is a column of the
  # design, so q <= degree + 1 <= T - 1. At the extreme q = T - 1 the two sides
  # are equal, never unequal:
  for (Tn in 3:8)
    expect_equal((Tn - 1) * Tn / 2 + Tn, Tn * (Tn + 1) / 2)

  # The check therefore guards the emission's contract rather than this
  # wrapper's, which matters because the same emission takes a non-polynomial
  # Lambda when factor mixture models arrive. The wrapper's own rule is what
  # rejects an over-parameterised trajectory.
  set.seed(109)
  Y <- matrix(rnorm(200 * 3), 200, 3)
  expect_error(fit_gmm(Y, times = 3, n_classes = 2, degree = 2),
               "would reproduce the occasion means exactly")

  # Just identified is allowed: 10 available, 6 for Psi and 4 residuals spent.
  Y4 <- matrix(rnorm(200 * 4), 200, 4)
  expect_no_error(
    suppressWarnings(
      fit_gmm(Y4, times = 4, n_classes = 2, degree = 2, random_effects = "all",
              residual = "occasion", n_init = 2, random_state = 1)))
})

test_that("estimation recovers a simulated growth mixture", {
  set.seed(2024)
  n <- 800; Tn <- 5
  cls   <- rbinom(n, 1, 0.5)
  icept <- ifelse(cls == 1, 4, 1) + rnorm(n, 0, 1.0)
  slope <- ifelse(cls == 1, 1.0, 0.1) + rnorm(n, 0, 0.3)
  Y <- outer(icept, rep(1, Tn)) + outer(slope, 0:(Tn - 1)) +
    matrix(rnorm(n * Tn, 0, 0.6), n, Tn)

  fit <- fit_gmm(Y, times = Tn, n_classes = 2, n_init = 10, random_state = 3)

  up <- which.max(fit$growth$means[, 2])
  expect_equal(unname(fit$growth$means[up, ]), c(4, 1.0), tolerance = 0.15)
  expect_equal(unname(fit$growth$means[-up, ]), c(1, 0.1), tolerance = 0.15)
  # Psi is class-equal by default, and the simulation used the same one.
  expect_equal(diag(fit$growth$psi[[1]]), c(1.0, 0.3^2), tolerance = 0.1,
               ignore_attr = TRUE)
  expect_equal(mean(fit$growth$residual_variance), 0.6^2, tolerance = 0.05)
  expect_true(fit$converged)
})

test_that("FIML keeps cases with unobserved occasions", {
  set.seed(113)
  n <- 800; Tn <- 5
  cls   <- rbinom(n, 1, 0.5)
  icept <- ifelse(cls == 1, 4, 1) + rnorm(n, 0, 1.0)
  slope <- ifelse(cls == 1, 1.0, 0.1) + rnorm(n, 0, 0.3)
  Y <- outer(icept, rep(1, Tn)) + outer(slope, 0:(Tn - 1)) +
    matrix(rnorm(n * Tn, 0, 0.6), n, Tn)
  Y[sample(seq_len(n), 150), 4] <- NA
  Y[sample(seq_len(n), 100), 5] <- NA

  fit <- fit_gmm(Y, times = Tn, n_classes = 2, n_init = 10, random_state = 3)

  expect_equal(nrow(fit$data), n)
  expect_equal(fit$missing_data$n_missing, 250L)
  up <- which.max(fit$growth$means[, 2])
  expect_equal(unname(fit$growth$means[up, ]), c(4, 1.0), tolerance = 0.2)
  expect_true(all(fit$growth$residual_variance > 0))
})

test_that("time_scores rescale the slope without moving the fit", {
  set.seed(127)
  Y <- .sim_gmm(n = 500)

  a <- fit_gmm(Y, times = 4, n_classes = 2, n_init = 5, random_state = 3)
  b <- fit_gmm(Y, times = 4, n_classes = 2, n_init = 5, random_state = 3,
               time_scores = c(0, 2, 4, 6))

  expect_equal(a$metrics$ll, b$metrics$ll, tolerance = 1e-4)
  expect_equal(a$growth$means[, 1], b$growth$means[, 1], tolerance = 1e-3)
  expect_equal(a$growth$means[, 2], 2 * b$growth$means[, 2], tolerance = 1e-3)
  # Halving the slope halves its standard deviation, so the slope variance goes
  # down by four and the covariance by two.
  expect_equal(a$growth$psi[[1]][2, 2], 4 * b$growth$psi[[1]][2, 2],
               tolerance = 1e-2)
  expect_equal(a$growth$psi[[1]][1, 2], 2 * b$growth$psi[[1]][1, 2],
               tolerance = 1e-2)
})

test_that("classes are sorted by size with all three parameter blocks together", {
  set.seed(131)
  Y <- .sim_gmm(n = 600, prob = 0.25)

  fit <- suppressWarnings(
    fit_gmm(Y, times = 4, n_classes = 2, psi = "free", residual_equal = FALSE,
            n_init = 10, random_state = 5))

  expect_true(fit$weights[1] >= fit$weights[2])
  # The larger class is the one at zero, so after sorting the first row of every
  # per-class block must belong to it. A block left unpermuted would pair the
  # big class's mean with the small class's variances.
  expect_lt(fit$growth$means[1, 1], fit$growth$means[2, 1])
  expect_lt(fit$growth$residual_variance[1, 1], 4)
  expect_equal(length(fit$growth$psi), 2L)
})

test_that("print and plot work, with and without random effects", {
  set.seed(137)
  Y <- .sim_gmm(n = 400)

  g <- fit_gmm(Y, times = 4, n_classes = 2, n_init = 5, random_state = 1)
  l <- fit_gmm(Y, times = 4, n_classes = 2, random_effects = "none",
               n_init = 5, random_state = 1)

  out_g <- paste(capture.output(print(g)), collapse = "\n")
  out_l <- paste(capture.output(print(l)), collapse = "\n")
  expect_match(out_g, "GROWTH FACTOR \\(CO\\)VARIANCE")
  expect_match(out_g, "held equal across classes")
  expect_match(out_g, "intercept, linear")
  expect_false(grepl("GROWTH FACTOR \\(CO\\)VARIANCE", out_l))
  expect_match(out_l, "none \\(latent class growth analysis\\)")

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot(g))
  expect_no_error(plot(g, observed = "cases"))
  expect_no_error(plot(g, observed = "none"))
})

test_that("the BLRT generator draws whole correlated trajectories", {
  set.seed(139)
  Y <- .sim_gmm(n = 800)
  fit <- fit_gmm(Y, times = 4, n_classes = 2, n_init = 5, random_state = 1)

  set.seed(11)
  N   <- 4000
  cl  <- rep(1L, N)
  gen <- generate_synthetic_data(fit$mm, cl, N)
  expect_equal(dim(gen), c(N, 4L))

  # The point of the generator for this emission: within a class the occasions
  # are correlated, because the random effects are what a case carries across
  # all of them. Generating occasion by occasion would give a diagonal
  # covariance and a null distribution for a model nobody fitted.
  expect_equal(cov(gen), .sn_sigma(fit$mm, 1L), tolerance = 0.15,
               ignore_attr = TRUE)
  expect_gt(cor(gen[, 1], gen[, 2]), 0.3)
  expect_equal(colMeans(gen), fit$growth$fitted[1, ], tolerance = 0.1,
               ignore_attr = TRUE)
})

test_that("predictors of class membership ride the existing machinery", {
  skip_on_cran()
  set.seed(31)
  Y <- .sim_gmm()
  x <- rnorm(nrow(Y))

  fit <- suppressMessages(
    fit_gmm(Y, times = 4, n_classes = 2, predictors = x, n_steps = 1,
            n_init = 10, random_state = 2026))

  # 12 measurement + 1 proportion becomes 12 + 2 regression parameters, since
  # the covariate model replaces the pooled class weights entirely.
  expect_equal(fit$metrics$n_params, 13L)
  expect_true(all(is.finite(fit$sm$parameters$beta)))
  expect_equal(unname(fit$sm$parameters$beta[2, ]), c(0, 0))
})
