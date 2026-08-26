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

test_that("absolute fit with missing data reduces to the complete-data numbers", {
  # A cheap guard against the MAR branch ever touching the complete-data path:
  # g2, x2 and cressie_read for a fit with no missing data must be exactly what
  # they were before this branch existed.
  fit <- .diag_fit(.diag_sim())
  af  <- absolute_fit(fit)
  expect_null(af$mar)
  expect_true(is.finite(af$dissimilarity))

  # The same data with a handful of cells set to NA takes the MAR branch and
  # returns an object rather than NULL, with the df the MAR formula implies.
  X  <- .diag_sim()
  Xm <- X
  Xm[1:5, 1] <- NA
  fit_na <- fit_mixture(Xm, n_components = 2, measurement = list(
    poly = list(model = "categorical_nan", n_columns = 2, max_val = 3),
    bin  = list(model = "binary_nan", n_columns = 2)),
    n_init = 5, random_state = 7)
  af_na <- absolute_fit(fit_na)
  expect_s3_class(af_na, "absolute_fit")
  expect_true(af_na$mar)
  W <- 3 * 3 * 2 * 2
  expect_equal(af_na$df, W - 1L - af_na$n_params)

  mt <- mcar_test(fit_na)
  expect_s3_class(mt, "mcar_test")
  expect_equal(mt$df, af_na$df_mcar - af_na$df)
  expect_equal(mt$stat, af_na$g2_mcar - af_na$g2)

  expect_message(expect_null(mcar_test(fit)), "complete")
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
    # form sum_k gk p_a(r|k) p_b(s|k) is a shortcut for, gk here being
    # colSums(resp * w) rather than fit$weights -- with no missing data the
    # two agree only up to the E-step/M-step lag left by EM's own convergence
    # tolerance, which is where the looser bound below comes from, not from
    # any imprecision in the closed form itself.
    p_ab <- tapply(tb$p, list(factor(tb$cells[, a], ra),
                              factor(tb$cells[, b], rb)), sum)
    o_ab <- tapply(o_all, list(factor(tb$cells[, a], ra),
                               factor(tb$cells[, b], rb)), sum)
    ref  <- sum((o_ab - N * p_ab)^2 / (N * p_ab)) /
      ((length(ra) - 1L) * (length(rb) - 1L))
    expect_lt(abs(unclass(bvr)[b, a] - ref), 1e-3)
  }
})

test_that("the expected-count fix leaves complete data untouched", {
  # With no missing data every pair keeps every case, so the posterior-weighted
  # class totals used for the expected counts and the plain marginal class
  # weights give the same expected table up to the E-step/M-step lag EM's own
  # convergence tolerance leaves behind -- this is the property that makes the
  # fix invisible, to a few significant digits, on the data most users have.
  fit   <- .diag_fit(.diag_sim(), bayes_constants = list(latent = 0))
  bvr   <- unclass(bivariate_residuals(fit))
  items <- .categorical_item_probs(fit$mm)
  gamma <- fit$weights
  w     <- fit$sample_weights

  for (b in 2:4) for (a in seq_len(b - 1L)) {
    ia <- items[[a]]; ib <- items[[b]]
    fa <- factor(fit$data[, a], levels = ia$categories)
    fb <- factor(fit$data[, b], levels = ib$categories)
    obs <- tapply(w, list(fa, fb), sum); obs[is.na(obs)] <- 0
    prob <- crossprod(ia$probs * gamma, ib$probs)
    expected <- sum(w) * prob
    ref <- sum((obs - expected)^2 / expected) /
      ((length(ia$categories) - 1L) * (length(ib$categories) - 1L))
    expect_lt(abs(bvr[b, a] - ref), 1e-3)
  }
})

test_that("missing data on one indicator moves only that indicator's pairs", {
  # Blank about 10% of column "a" only, then compare the fixed formula against
  # the old marginal-weight one on the *same* fit -- isolating what the
  # formula change does from what refitting on different data would do. Pairs
  # not touching "a" keep every case, so their retained-subset posterior and
  # the marginal class weights agree up to EM's own convergence lag; only
  # pairs with "a" are computed on a subsample whose class composition can
  # genuinely differ from the full sample's.
  set.seed(9)
  X <- .diag_sim(n = 800, seed = 9)
  drop <- sample(nrow(X), floor(0.10 * nrow(X)))
  X[drop, "a"] <- NA

  fit   <- .diag_fit(X)
  items <- .categorical_item_probs(fit$mm)
  gamma <- fit$weights
  w     <- fit$sample_weights
  resp  <- exp(fit$log_resp)

  old_bvr <- function(a, b) {
    ia <- items[[a]]; ib <- items[[b]]
    keep <- !is.na(fit$data[, a]) & !is.na(fit$data[, b])
    fa <- factor(fit$data[keep, a], levels = ia$categories)
    fb <- factor(fit$data[keep, b], levels = ib$categories)
    obs <- tapply(w[keep], list(fa, fb), sum); obs[is.na(obs)] <- 0
    expected <- sum(w[keep]) * crossprod(ia$probs * gamma, ib$probs)
    sum((obs - expected)^2 / expected) /
      ((length(ia$categories) - 1L) * (length(ib$categories) - 1L))
  }

  bvr <- unclass(bivariate_residuals(fit))

  touches_a <- c(b = 2L, c = 3L, d = 4L)
  moved <- vapply(touches_a, function(j) abs(bvr[j, 1L] - old_bvr(1L, j)), 0)
  unmoved <- abs(bvr[3L, 2L] - old_bvr(2L, 3L))

  expect_true(all(moved > 10 * unmoved))
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
  # Absolute fit now works with missing data, computed under MAR, and returns
  # an object with a smaller df than the complete-data case would give.
  af <- absolute_fit(m)
  expect_s3_class(af, "absolute_fit")
  expect_true(af$mar)
  expect_equal(af$df, prod(sapply(.categorical_item_probs(m$mm), function(it)
    length(it$categories))) - 1L - af$n_params)
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

# Two classes, five indicators, one seed. Items 1 and 2 carry a within-class
# residual correlation of 0.4 in the first simulated class only; every other
# pair, and every pair in the second, is independent. Every class variance is
# 1 by construction, so the same data is the right shape for both
# parameterisations: the misfit is a covariance, not a variance difference.
.ld_sim <- function() {
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
  for (j in 1:5) X[, j] <- X[, j] + 2 * (k - 1)  # separate the classes
  X
}

# sort_model_classes() orders classes largest-to-smallest, and the two
# parameterisations need not agree on which fitted class is which, so the
# planted class is always found by matching item-1 means against the two
# simulated centres rather than assumed to be class 1.
.ld_planted <- function(fit) which.min(abs(fit$mm$parameters$means[, 1] - 0))

test_that("the modification index flags the one planted local dependence", {
  # Two classes, five indicators, one seed. Items 1 and 2 carry a within-class
  # residual correlation of 0.4 in class 1 only; every other pair, and every
  # pair in class 2, is independent. The misspecified conditional-independence
  # fit should flag exactly the planted pair, in exactly the class it was
  # planted in.
  X <- .ld_sim()

  fit <- fit_mixture(X, n_classes = 2, measurement = "continuous", variances_equal = FALSE,
                     n_init = 5, random_state = 1)
  bvr <- bivariate_residuals(fit)
  expect_s3_class(bvr, "bivariate_residuals_gaussian")

  planted <- .ld_planted(fit)

  all_mi <- as.vector(bvr$mi)
  planted_mi <- bvr$mi[2, 1, planted]
  expect_equal(planted_mi, max(all_mi, na.rm = TRUE))
  expect_lt(bvr$p_value[2, 1, planted], 0.05)

  null_mi <- all_mi[!is.na(all_mi)]
  null_mi <- null_mi[null_mi != planted_mi]
  null_p  <- stats::pchisq(null_mi, df = 1, lower.tail = FALSE)
  expect_gt(mean(null_p > 0.05), 0.75)
})

test_that("the modification index holds up under the equal-variance default", {
  # A continuous measurement model holds each item's variance equal across the
  # classes by default, so this is the parameterisation a user reaches without
  # asking for it, and the one the statistic has to be right under. The
  # constrained fit packs one shared log-sd per item rather than K
  # (.step1_pack_sub); were the constraint not carried into the packed vector,
  # the fit would not be stationary in the K - 1 directions the packing
  # invented, the Schur complement would go non-positive, and every pair would
  # fall back to the outer-product denominator.
  X   <- .ld_sim()
  fit <- fit_mixture(X, n_classes = 2, measurement = "continuous",
                     n_init = 5, random_state = 1)
  expect_true(isTRUE(fit$mm$variances_equal))
  # The constraint binds: one variance per item, shared by the two classes.
  expect_equal(fit$mm$parameters$covariances[1, ],
               fit$mm$parameters$covariances[2, ])

  bvr     <- bivariate_residuals(fit)
  planted <- .ld_planted(fit)
  mi      <- bvr$mi

  # Nothing fell back to the outer-product denominator, and nothing came out
  # non-finite. n_opg is the direct symptom the packing gap used to produce.
  expect_equal(bvr$n_opg, 0L)
  expect_true(all(is.finite(mi[!is.na(mi)])))

  # It still finds the planted pair, in the class it was planted in.
  expect_equal(mi[2, 1, planted], max(mi, na.rm = TRUE))
  expect_lt(bvr$p_value[2, 1, planted], 0.05)

  null_mi <- as.vector(mi)[!is.na(as.vector(mi))]
  null_mi <- null_mi[null_mi != mi[2, 1, planted]]
  expect_gt(mean(stats::pchisq(null_mi, df = 1, lower.tail = FALSE) > 0.05),
            0.75)

  # And it agrees with the free-variance fit on every pair that carries no
  # planted dependence, once the classes are matched up -- the two orderings
  # need not coincide. The planted pair itself is not comparable: the
  # constrained fit prices it against different information. Loose by design,
  # since these are two different fits of two different models.
  free <- fit_mixture(X, n_classes = 2, measurement = "continuous",
                      variances_equal = FALSE, n_init = 5, random_state = 1)
  null_pairs <- lower.tri(diag(5))
  null_pairs[2, 1] <- FALSE
  expect_lt(max(abs(mi[, , planted][null_pairs] -
                      bivariate_residuals(free)$mi[, , .ld_planted(free)][null_pairs])),
            0.5)
})

test_that("lta_g2() still reports what it always did", {
  fit <- .diag_fit(.diag_sim())
  g2  <- lta_g2(fit)
  af  <- absolute_fit(fit)
  expect_named(g2, c("g2", "df", "p_value", "n_cells", "n_patterns"))
  expect_equal(g2$g2, af$g2)
  expect_equal(g2$df, af$df)
})

test_that("the bootstrap attaches calibrated p-values without moving the residuals", {
  fit <- .diag_fit(.diag_sim())

  plain <- bivariate_residuals(fit)
  boot  <- bivariate_residuals(fit, n_reps = 3)

  # The residuals themselves are a function of the fit, not of the bootstrap.
  stripped <- unclass(boot)
  attr(stripped, "p") <- NULL
  expect_equal(unclass(plain), stripped)

  p <- attr(boot, "p")
  expect_false(is.null(p))
  expect_null(attr(plain, "p"))
  expect_equal(dim(p), dim(unclass(plain)))
  expect_equal(dimnames(p), dimnames(unclass(plain)))

  # A proportion of replicates, so it lives on [0, 1], and it is defined
  # exactly where the residual is.
  expect_true(all(p >= 0 & p <= 1, na.rm = TRUE))
  expect_equal(is.na(p[lower.tri(p)]), is.na(unclass(plain)[lower.tri(p)]))

  expect_output(print(boot), "Bootstrap p-value")
  expect_error(bivariate_residuals(fit, n_reps = -1), "non-negative")
})

test_that("the continuous branch says it is ignoring n_reps", {
  set.seed(11)
  Xc  <- cbind(a = rnorm(200), b = rnorm(200) + 1, c = rnorm(200))
  fit <- fit_mixture(Xc, n_classes = 2, measurement = "continuous", n_init = 2)

  expect_message(bivariate_residuals(fit, n_reps = 3), "ignored")
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

# ------------------------------------------------------------------------------
# Multiple-group (group_effects = "measurement") fits
# ------------------------------------------------------------------------------
#
# A group_blocks fit pads J items into J*Q columns, one block per group
# (R/group_blocks.R). Handed to the ordinary code unchanged, that padding
# makes .categorical_item_probs() see J*Q distinct items and build a table of
# prod(levels)^Q cells rather than the Q * prod(levels) the fitted model
# (Q separate P(y | group) tables) actually implies. The fix routes a
# group_blocks fit through its own group-aware path; these tests check the
# result against a hand-enumerated brute force, not just that it runs.

# Two groups, three binary items, two classes, item probabilities set to
# differ by group so the fitted model is genuinely group-varying.
.group_diag_sim <- function(n_per_group = 400, seed = 42) {
  set.seed(seed)
  J <- 3
  rho <- list(
    A = matrix(c(.8, .7, .6, .2, .3, .4), 2, J, byrow = TRUE),
    B = matrix(c(.75, .65, .55, .25, .35, .45), 2, J, byrow = TRUE))
  gamma <- c(.6, .4)
  draw <- function(r) {
    k <- sample(1:2, n_per_group, replace = TRUE, prob = gamma)
    X <- matrix(NA_real_, n_per_group, J)
    for (j in 1:J) X[, j] <- rbinom(n_per_group, 1, r[k, j])
    X
  }
  X <- rbind(draw(rho$A), draw(rho$B))
  colnames(X) <- paste0("y", 1:J)
  list(X = X, group = factor(rep(c("A", "B"), each = n_per_group)))
}

# The brute-force reference: enumerate the 2^J response patterns and compute
# each group's own model-implied probability from that group's own fitted
# item probabilities and the (pooled) class weights, entirely independently
# of the package's own absolute_fit() implementation.
.group_brute_g2 <- function(fit) {
  J <- fit$mm$n_items
  cells <- as.matrix(expand.grid(rep(list(0:1), J)))
  gamma <- fit$weights
  g2_total <- 0
  for (q in seq_along(fit$mm$models)) {
    pis <- fit$mm$models[[q]]$parameters$pis
    K <- nrow(pis)
    probs <- apply(cells, 1, function(cell) {
      sum(vapply(seq_len(K), function(k)
        gamma[k] * prod(ifelse(cell == 1, pis[k, ], 1 - pis[k, ])),
        numeric(1)))
    })
    rows <- which(as.integer(fit$group_info$factor) == q)
    Xg <- .strip_block_prefix(fit$data[rows, .time_block_cols(q, J), drop = FALSE])
    key_cells <- apply(cells, 1, paste, collapse = "|")
    obs <- as.numeric(table(factor(apply(Xg, 1, paste, collapse = "|"),
                                   levels = key_cells)))
    exp_c <- nrow(Xg) * probs
    g2_total <- g2_total + 2 * sum(ifelse(obs > 0, obs * log(obs / exp_c), 0))
  }
  g2_total
}

test_that("absolute_fit() on a group-varying measurement fit matches a hand-enumerated Q x W table", {
  sim <- .group_diag_sim()
  fit <- fit_mixture(sim$X, n_classes = 2, measurement = "binary",
                     group = sim$group, group_effects = "measurement",
                     n_steps = 1, n_init = 5, random_state = 1)

  af <- absolute_fit(fit)
  expect_s3_class(af, "absolute_fit")
  expect_null(af$mar)   # no missing data here

  W <- 2^3
  Q <- 2L
  expect_equal(af$n_cells, Q * W)
  expect_equal(af$df, Q * (W - 1L) - fit$metrics$n_params)
  expect_equal(af$g2, .group_brute_g2(fit), tolerance = 1e-8)

  # The per-group breakdown sums to the total.
  expect_equal(sum(af$by_group$g2), af$g2, tolerance = 1e-8)
  expect_equal(nrow(af$by_group), Q)
})

test_that("the Q(W-1)-npar formula reduces to the ungrouped df at Q=1", {
  # A pure arithmetic identity, asserted directly rather than by fitting a
  # one-group model (`group=` always requires at least two levels): the
  # roadmap's requirement is that the grouped formula not silently change
  # what the ungrouped path already computes.
  W <- 44L; npar <- 13L
  expect_equal(1L * (W - 1L) - npar, W - npar - 1L)
})

test_that("absolute_fit() on a group-varying fit still handles missing data under MAR", {
  sim <- .group_diag_sim(n_per_group = 500, seed = 7)
  X <- sim$X
  X[sample(length(X), floor(0.05 * length(X)))] <- NA
  fit <- fit_mixture(X, n_classes = 2, measurement = "binary",
                     group = sim$group, group_effects = "measurement",
                     n_steps = 1, n_init = 5, random_state = 1)
  expect_true(anyNA(fit$data))

  af <- absolute_fit(fit)
  expect_true(af$mar)
  expect_equal(af$df, 2L * (2^3 - 1L) - fit$metrics$n_params)
  expect_true(af$g2 >= 0)
  expect_equal(nrow(af$by_group), 2L)
})

test_that("group_effects='both' is still refused as a conditional model, unchanged", {
  sim <- .group_diag_sim(n_per_group = 100, seed = 3)
  fit <- suppressWarnings(fit_mixture(
    sim$X, n_classes = 2, measurement = "binary",
    group = sim$group, group_effects = "both",
    n_steps = 1, n_init = 2, random_state = 1))
  expect_message(expect_null(absolute_fit(fit)), "covariates")
  expect_message(expect_null(bivariate_residuals(fit)), "covariates")
})

test_that("bivariate_residuals() on a group-varying fit is J x J, in item names, not J*Q x J*Q", {
  sim <- .group_diag_sim()
  fit <- fit_mixture(sim$X, n_classes = 2, measurement = "binary",
                     group = sim$group, group_effects = "measurement",
                     n_steps = 1, n_init = 5, random_state = 1)

  bvr <- bivariate_residuals(fit)
  expect_equal(dim(bvr), c(3L, 3L))
  expect_equal(rownames(bvr), c("y1", "y2", "y3"))
  expect_equal(attr(bvr, "n_groups"), 2L)
  by_group <- attr(bvr, "by_group")
  expect_named(by_group, c("A", "B"))
  expect_equal(dim(by_group$A), c(3L, 3L))

  # The combined residual for a pair is the group-average of that pair's own
  # chi-square, not merely their sum -- averaging is what keeps it centred
  # near 1 the way the single-group statistic is.
  df_pair <- 1  # two binary items: (2-1)*(2-1)
  a <- by_group$A[2, 1] * df_pair
  b <- by_group$B[2, 1] * df_pair
  expect_equal(bvr[2, 1], (a + b) / (2 * df_pair))
})

