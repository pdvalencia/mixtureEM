# The one-score class: the package's single definition of a degenerate fit
# (R/one_score.R). What it tests, per family, is the variances a class owns;
# what it asks of each is whether the class sits on one value of the data that
# variance governs.

test_that("only a variance a class owns is tested, pooled over what it governs", {
  K <- 2; J <- 3
  free <- gaussian_model(K, "gaussian_diag")
  free$parameters$means <- matrix(0, K, J)
  expect_length(.dispersion_groups(free), K * J)

  equal <- gaussian_model(K, "gaussian_diag", variances_equal = TRUE)
  equal$parameters$means <- matrix(0, K, J)
  expect_length(.dispersion_groups(equal), 0L)

  # Item 1 held equal over three occasions: one variance per class for it,
  # governing all three of its columns; items 2 and 3 stay one per occasion.
  tb <- time_blocks_model(K, J, 3, sub_model = "gaussian_diag",
                          invariant_items = 1L)
  for (b in 1:3) tb$models[[b]]$parameters$means <- matrix(0, K, J)
  gs <- .dispersion_groups(tb)
  expect_length(gs, K * (1 + 2 * 3))
  pooled <- Filter(function(g) length(g$cols) > 1L, gs)
  expect_length(pooled, K)
  expect_equal(sort(pooled[[1]]$cols), c(1L, 4L, 7L))
  expect_equal(sort(pooled[[1]]$blocks), 1:3)

  L <- .lcga_design(0:3, 1)
  lc <- lcga_model(K, L, family = "gaussian")
  expect_length(.dispersion_groups(lc), K)
  expect_equal(.dispersion_groups(lc)[[1]]$cols, 1:4)
  expect_length(.dispersion_groups(lcga_model(K, L, family = "binomial")), 0L)

  expect_length(.dispersion_groups(structured_normal_model(
    K, L, psi = "equal", residual = "occasion", residual_equal = TRUE)), 0L)
  expect_length(.dispersion_groups(structured_normal_model(
    K, L, psi = "free", residual = "occasion", residual_equal = TRUE)), K * 4)
  expect_length(.dispersion_groups(structured_normal_model(
    K, L, psi = "equal", residual = "constant", residual_equal = FALSE)), K)
})

test_that("a class is flagged only when it sits on one value of all its variance governs", {
  K <- 2; n <- 100
  tb <- time_blocks_model(K, 1, 2, sub_model = "gaussian_diag",
                          invariant_items = 1L)
  for (b in 1:2) tb$models[[b]]$parameters$means <- matrix(0, K, 1)
  set.seed(1)
  X <- cbind(c(rep(5, 30), rnorm(70)), rnorm(n))
  # Class 1 is the thirty tied cases at occasion 1 -- but the same variance
  # also governs occasion 2, where those cases are spread out.
  on_ties <- c(rep(1, 30), rep(0, 70))
  post <- list(cbind(on_ties, 1 - on_ties), cbind(on_ties, 1 - on_ties))
  expect_null(.one_score_cells(tb, X, post))

  # Free by occasion, the occasion-1 variance is the class's alone.
  tb_free <- time_blocks_model(K, 1, 2, sub_model = "gaussian_diag")
  for (b in 1:2) tb_free$models[[b]]$parameters$means <- matrix(0, K, 1)
  hit <- .one_score_cells(tb_free, X, post)
  expect_equal(nrow(hit), 1L)
  expect_equal(hit$class, 1L); expect_equal(hit$col, 1L); expect_equal(hit$value, 5)

  # Per-occasion posteriors: at occasion 2 the status holding the tied cases
  # at occasion 1 is someone else, so nothing at occasion 2 is flagged.
  post2 <- list(cbind(on_ties, 1 - on_ties), cbind(rep(0.5, n), rep(0.5, n)))
  expect_equal(.one_score_cells(tb_free, X, post2)$col, 1L)

  # Case weights count: weighting the tied cases down below 90% of the class
  # clears the flag.
  mixed <- c(rep(1, 30), rep(1, 2), rep(0, 68))     # 30 of 32 tied: 94%
  w <- c(rep(0.1, 30), rep(1, 70))
  free <- gaussian_model(K, "gaussian_diag")
  free$parameters$means <- matrix(0, K, 2)
  expect_false(is.null(.one_score_cells(free, X, cbind(mixed, 1 - mixed))))
  expect_null(.one_score_cells(free, X, cbind(mixed, 1 - mixed), weights = w))
})

test_that("a latent profile search sets aside a start whose class sits on a ceiling", {
  skip_on_cran()
  # A strong ceiling on a rounded scale, three true profiles fitted with four
  # and no prior: some starts let a class sit on the ceiling and run off to an
  # unbounded likelihood. They are set aside; the fit reports an admissible
  # solution and no collapse.
  set.seed(6181)
  n <- 600
  z <- sample(1:3, n, TRUE, c(.45, .35, .2))
  mu <- rbind(c(2.4, 2.6, 2.2), c(3.4, 3.6, 3.3), c(4.9, 4.9, 4.3))
  y <- mu[z, ] + matrix(rnorm(3 * n), n) * c(.55, .45, .35)[z]
  y <- pmin(pmax(round(y * 4) / 4, 1), 5)
  colnames(y) <- c("a", "b", "c")
  off <- list(categorical = 0, latent = 0, poisson = 0, variances = 0)

  fit <- suppressWarnings(fit_mixture(y, n_classes = 4, measurement = "continuous",
    variances_equal = FALSE, n_init = 20, random_state = 1,
    bayes_constants = off))
  expect_gt(fit$metrics$n_set_aside, 0L)
  expect_null(fit$degenerate)
  expect_null(.one_score_cells(fit$mm, y, exp(fit$log_resp)))
  expect_lt(fit$metrics$ll, max(fit$start_lls))

  # Equal variances bound the likelihood: nothing to set aside.
  eq <- suppressWarnings(fit_mixture(y, n_classes = 4, measurement = "continuous",
    n_init = 5, random_state = 1, bayes_constants = off))
  expect_equal(eq$metrics$n_set_aside %||% 0L, 0L)
})
