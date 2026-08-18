# ==============================================================================
# Absolute fit, bivariate residuals and the classification table
# ==============================================================================
#
# Two kinds of check, and the first kind matters more than it looks.
#
# The implementations deliberately avoid enumerating the response-pattern
# table: the Pearson X^2 folds every unobserved cell into a single closing
# term, and the bivariate residuals take the two-way margin from conditional
# independence in closed form. Both are algebraic shortcuts, and a shortcut
# that is subtly wrong still returns a plausible number. So the first block
# below re-derives each statistic the slow, obvious way on a table small
# enough to enumerate, and requires the two to agree to machine precision.
# These need no external data and run everywhere.
#
# Checks that anchor the same statistics against external reference output
# live in the internal validation suite, which is not part of this package.

# ------------------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------------------

# Four categorical indicators: two three-category, two binary, matching the
# shape of that example so the same fixtures serve both blocks.
.diag_sim <- function(n = 400, seed = 42) {
  set.seed(seed)
  k <- sample(1:2, n, replace = TRUE, prob = c(0.6, 0.4))
  p3 <- list(matrix(c(.7, .2, .1, .1, .3, .6), 2, 3, byrow = TRUE),
             matrix(c(.2, .5, .3, .6, .2, .2), 2, 3, byrow = TRUE))
  draw3 <- function(pr) vapply(k, function(kk) sample(1:3, 1, prob = pr[kk, ]),
                               integer(1))
  cbind(a = draw3(p3[[1]]), b = draw3(p3[[2]]),
        c = rbinom(n, 1, c(.8, .3)[k]), d = rbinom(n, 1, c(.25, .7)[k]))
}

.diag_spec <- list(poly = list(model = "categorical", n_columns = 2, max_val = 3),
                   bin  = list(model = "binary",      n_columns = 2))

.diag_fit <- function(X, K = 2, ...)
  fit_mixture(X, n_components = K, measurement = .diag_spec,
              n_init = 5, random_state = 7, ...)

# Every cell of the response-pattern table, in data-column order, with the
# probability the model gives it. The slow reference the shortcuts are checked
# against.
.diag_full_table <- function(fit) {
  items <- .categorical_item_probs(fit$mm)
  grid  <- as.matrix(expand.grid(lapply(items, `[[`, "categories")))
  colnames(grid) <- colnames(fit$data)
  ll <- log_likelihood(fit$mm, grid)
  list(cells = grid, p = as.vector(exp(ll) %*% fit$weights))
}

# Weighted observed count of each enumerated cell.
.diag_observed <- function(fit, cells) {
  w   <- fit$sample_weights
  key <- function(M) apply(M, 1, paste, collapse = "\r")
  obs <- tapply(w, factor(key(fit$data), levels = key(cells)), sum)
  ifelse(is.na(obs), 0, obs)
}

# ------------------------------------------------------------------------------
# The shortcuts against brute force
# ------------------------------------------------------------------------------

test_that("the model-implied table is a probability distribution", {
  fit <- .diag_fit(.diag_sim())
  tb  <- .diag_full_table(fit)
  expect_equal(sum(tb$p), 1, tolerance = 1e-12)
  expect_equal(nrow(tb$cells), 3 * 3 * 2 * 2)
})

test_that("X2, G2 and Cressie-Read match their enumerated definitions", {
  fit <- .diag_fit(.diag_sim())
  af  <- absolute_fit(fit)
  tb  <- .diag_full_table(fit)

  N <- sum(fit$sample_weights)
  o <- .diag_observed(fit, tb$cells)
  e <- N * tb$p

  # Pearson over every cell, including the ones no case took. This is the term
  # absolute_fit() replaces with (N - sum of observed expected counts).
  expect_equal(af$x2, sum((o - e)^2 / e), tolerance = 1e-10)

  # G2 and Cressie-Read weight each cell by its observed count, so unobserved
  # cells drop out; enumerating them must change nothing.
  pos <- o > 0
  expect_equal(af$g2, 2 * sum(o[pos] * log(o[pos] / e[pos])), tolerance = 1e-10)
  lambda <- 2 / 3
  expect_equal(af$cressie_read,
               2 / (lambda * (lambda + 1)) *
                 sum(o[pos] * ((o[pos] / e[pos])^lambda - 1)),
               tolerance = 1e-10)

  # Cressie-Read sits between the other two, since lambda = 2/3 lies between
  # G2's limit at 0 and Pearson's 1. Which of the two it is closer to, and even
  # which way round they fall, depends on the table: on these data the family
  # decreases in lambda (20.211, 19.451, 19.235), on the benchmark data used in
  # the internal validation suite
  # it increases (22.087, 22.662, 23.511). Only the bracketing is a property of
  # the statistic rather than of the data.
  expect_gt(af$cressie_read, min(af$g2, af$x2))
  expect_lt(af$cressie_read, max(af$g2, af$x2))

  expect_equal(af$n_cells, 36L)
  expect_equal(af$df, 36L - af$n_params - 1L)
})

test_that("the closed-form bivariate margin equals the enumerated one", {
  fit   <- .diag_fit(.diag_sim())
  bvr   <- bivariate_residuals(fit)
  tb    <- .diag_full_table(fit)
  N     <- sum(fit$sample_weights)
  o_all <- .diag_observed(fit, tb$cells)

  items <- .categorical_item_probs(fit$mm)
  for (b in 2:4) for (a in seq_len(b - 1L)) {
    ra <- items[[a]]$categories; rb <- items[[b]]$categories
    # Marginalise the full table down to this pair, the definition the closed
    # form sum_k gamma_k p_a(r|k) p_b(s|k) is a shortcut for.
    p_ab <- tapply(tb$p, list(factor(tb$cells[, a], ra),
                              factor(tb$cells[, b], rb)), sum)
    o_ab <- tapply(o_all, list(factor(tb$cells[, a], ra),
                               factor(tb$cells[, b], rb)), sum)
    ref  <- sum((o_ab - N * p_ab)^2 / (N * p_ab)) /
      ((length(ra) - 1L) * (length(rb) - 1L))
    expect_equal(unclass(bvr)[b, a], ref, tolerance = 1e-10,
                 label = sprintf("BVR[%d,%d]", b, a))
  }
})

test_that("a one-class model's BVR is the ordinary test of independence", {
  # With one class the model-implied joint is the product of the fitted
  # marginals, and for K = 1 the Dirichlet prior is exactly neutral -- it adds
  # one pseudo-observation at the marginal it is already estimating -- so the
  # fitted marginals are the observed ones and the residual reduces to
  # Pearson's chi-square for a two-way table. An independent closed-form
  # target with no rounding in it anywhere.
  X   <- .diag_sim()
  fit <- .diag_fit(X, K = 1)
  bvr <- unclass(bivariate_residuals(fit))

  for (b in 2:4) for (a in seq_len(b - 1L)) {
    tt  <- table(X[, a], X[, b])
    ref <- unname(stats::chisq.test(tt, correct = FALSE)$statistic) /
      ((nrow(tt) - 1L) * (ncol(tt) - 1L))
    expect_equal(bvr[b, a], ref, tolerance = 1e-9,
                 label = sprintf("BVR[%d,%d]", b, a))
  }
})

test_that("bivariate residuals honour case weights", {
  # Collapsing duplicated rows into frequency weights must not move anything:
  # the same data, told two ways.
  X    <- .diag_sim(n = 600, seed = 3)
  key  <- apply(X, 1, paste, collapse = "\r")
  Xw   <- X[!duplicated(key), , drop = FALSE]
  freq <- as.vector(table(factor(key, levels = key[!duplicated(key)])))

  long  <- .diag_fit(X)
  short <- fit_mixture(Xw, n_components = 2, measurement = .diag_spec,
                       weights = freq, weight_type = "frequency",
                       n_init = 5, random_state = 7)

  expect_equal(long$metrics$ll, short$metrics$ll, tolerance = 1e-6)
  expect_equal(unclass(bivariate_residuals(long)),
               unclass(bivariate_residuals(short)), tolerance = 1e-6)
  expect_equal(absolute_fit(long)$g2, absolute_fit(short)$g2, tolerance = 1e-6)
  expect_equal(attr(classification_table(long), "error"),
               attr(classification_table(short), "error"), tolerance = 1e-6)
})

test_that("the classification table reconciles both sets of margins", {
  fit <- .diag_fit(.diag_sim())
  ct  <- classification_table(fit)
  tab <- unclass(ct)

  resp  <- exp(fit$log_resp)
  w     <- fit$sample_weights
  modal <- max.col(resp, ties.method = "first")

  # Columns are the modal counts, rows the model-expected class sizes; the two
  # disagree by exactly what modal assignment does to the class proportions.
  expect_equal(unname(colSums(tab)),
               as.vector(tapply(w, factor(modal, 1:2), sum)), tolerance = 1e-10)
  expect_equal(unname(rowSums(tab)), as.vector(colSums(resp * w)),
               tolerance = 1e-10)
  expect_equal(sum(tab), sum(w), tolerance = 1e-10)
  expect_equal(attr(ct, "error"), 1 - sum(diag(tab)) / sum(w))

  # A perfectly separated model classifies without error.
  set.seed(11)
  Xsep <- rbind(matrix(rep(c(1, 1, 0, 0), each = 100), 100, 4),
                matrix(rep(c(3, 3, 1, 1), each = 100), 100, 4))
  colnames(Xsep) <- c("a", "b", "c", "d")
  expect_lt(attr(classification_table(.diag_fit(Xsep)), "error"), 1e-6)
})

# ------------------------------------------------------------------------------
# Where the statistics refuse to apply
# ------------------------------------------------------------------------------

test_that("the diagnostics decline models they are not defined for", {
  set.seed(5)
  Y <- matrix(rnorm(400), ncol = 4)
  g <- fit_mixture(Y, n_components = 2, measurement = "continuous",
                   variances_equal = FALSE, n_init = 2)
  expect_message(expect_null(absolute_fit(g)), "categorical")
  # A plain continuous measurement model gets the modification-index
  # statistic instead of a refusal.
  expect_s3_class(bivariate_residuals(g), "bivariate_residuals_gaussian")
  # The classification table reads only the posterior, so it applies here.
  expect_s3_class(classification_table(g), "classification_table")

  # A continuous model with missing data has no bivariate-normal augmented
  # density defined here, and neither does one nested inside a mixed or
  # repeated-measures measurement model.
  Ym <- Y; Ym[1, 1] <- NA
  gm <- fit_mixture(Ym, n_components = 2,
                    measurement = "continuous_nan", n_init = 2)
  expect_message(expect_null(bivariate_residuals(gm)), "missing data")

  X <- .diag_sim(n = 200)
  Xm <- X; Xm[1:10, 2] <- NA
  m <- fit_mixture(Xm, n_components = 2, n_init = 2, random_state = 1,
                   measurement = list(
                     poly = list(model = "categorical_nan", n_columns = 2,
                                 max_val = 3),
                     bin  = list(model = "binary_nan", n_columns = 2)))
  expect_message(expect_null(absolute_fit(m)), "missing data")
  # Bivariate residuals survive missingness: each pair uses the cases that
  # observe both items.
  expect_false(anyNA(unclass(bivariate_residuals(m))[lower.tri(diag(4))]))

  # Covariates make the implied table case-specific.
  set.seed(6)
  Z <- cbind(1, rnorm(nrow(X)))
  cv <- fit_mixture(X, Y = Z, n_components = 2, measurement = .diag_spec,
                    structural = "covariate", n_steps = 1, n_init = 2)
  expect_message(expect_null(absolute_fit(cv)), "covariates")
  expect_message(expect_null(bivariate_residuals(cv)), "covariates")
})

# ------------------------------------------------------------------------------
# Local dependence, continuous indicators
# ------------------------------------------------------------------------------

test_that("the modification index flags the one planted local dependence", {
  # Two classes, five indicators, one seed. Items 1 and 2 carry a within-class
  # residual correlation of 0.4 in class 1 only; every other pair, and every
  # pair in class 2, is independent. The misspecified conditional-independence
  # fit should flag exactly the planted pair, in exactly the class it was
  # planted in.
  set.seed(123)
  n <- 800
  k <- sample(1:2, n, replace = TRUE)
  X <- matrix(rnorm(n * 5), n, 5)
  # Induce the class-1 residual correlation between items 1 and 2 by mixing in
  # a shared latent shock, scaled so cor(X[,1], X[,2]) = 0.4 among class 1.
  shock <- rnorm(n)
  rho <- 0.4
  in1 <- k == 1
  X[in1, 1] <- sqrt(rho) * shock[in1] + sqrt(1 - rho) * X[in1, 1]
  X[in1, 2] <- sqrt(rho) * shock[in1] + sqrt(1 - rho) * X[in1, 2]
  X[, 1] <- X[, 1] + 2 * (k - 1)   # separate the two classes on every item
  X[, 2] <- X[, 2] + 2 * (k - 1)
  X[, 3] <- X[, 3] + 2 * (k - 1)
  X[, 4] <- X[, 4] + 2 * (k - 1)
  X[, 5] <- X[, 5] + 2 * (k - 1)

  fit <- fit_mixture(X, n_components = 2, measurement = "continuous", variances_equal = FALSE,
                     n_init = 5, random_state = 1)
  bvr <- bivariate_residuals(fit)
  expect_s3_class(bvr, "bivariate_residuals_gaussian")

  # sort_model_classes() orders classes largest-to-smallest, so find which
  # fitted class corresponds to the planted (originally "class 1") group by
  # matching item-1 means against the two simulated centres.
  planted <- which.min(abs(fit$mm$parameters$means[, 1] - 0))

  all_mi <- as.vector(bvr$mi)
  planted_mi <- bvr$mi[2, 1, planted]
  expect_equal(planted_mi, max(all_mi, na.rm = TRUE))
  expect_lt(bvr$p_value[2, 1, planted], 0.05)

  null_mi <- all_mi[!is.na(all_mi)]
  null_mi <- null_mi[null_mi != planted_mi]
  null_p  <- stats::pchisq(null_mi, df = 1, lower.tail = FALSE)
  expect_gt(mean(null_p > 0.05), 0.75)
})

test_that("lta_g2() still reports what it always did", {
  fit <- .diag_fit(.diag_sim())
  g2  <- lta_g2(fit)
  af  <- absolute_fit(fit)
  expect_named(g2, c("g2", "df", "p_value", "n_cells", "n_patterns"))
  expect_equal(g2$g2, af$g2)
  expect_equal(g2$df, af$df)
})

test_that("the diagnostics print without error", {
  fit <- .diag_fit(.diag_sim())
  expect_output(print(absolute_fit(fit)), "L-squared")
  expect_output(print(bivariate_residuals(fit)), "Largest")
  expect_output(print(classification_table(fit)), "Classification error")
  expect_output(cd <- classification_diagnostics(fit), "AvePP")
  expect_named(cd, c("ave_pp", "table", "error"))

  # The heat table, including a residual past the scale's saturation point:
  # clamping, not zlim, is what has to keep that cell drawn.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(bivariate_residuals(fit)))
  big <- bivariate_residuals(fit)
  big[4, 2] <- 99
  expect_silent(plot(big, max_shade = 4))
})

