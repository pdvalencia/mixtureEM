# A one-step fit with covariates on class membership used to refuse the
# response-pattern economy outright. It no longer does: the key widens to
# (response pattern, covariate values): two cases are interchangeable only
# when responses and covariates both agree. The claim is again a saving and
# not a change, and these tests are what stands behind it.

.cc_sim <- function(seed = 11, n = 1500) {
  set.seed(seed)
  # Two discrete covariates, so cases really do share covariate values. The
  # class a case belongs to depends on them; the responses depend on the class.
  z1 <- stats::rbinom(n, 1L, 0.5)
  z2 <- stats::rbinom(n, 2L, 0.5)
  eta <- -0.6 + 1.2 * z1 + 0.5 * z2
  s   <- 1L + stats::rbinom(n, 1L, 1 / (1 + exp(-eta)))
  rho <- matrix(c(.85, .80, .20, .25, .15,
                  .20, .25, .80, .75, .85), 2, 5, byrow = TRUE)
  X <- matrix(stats::rbinom(n * 5, 1L, rho[s, ]), n, 5)
  colnames(X) <- paste0("i", 1:5)
  list(X = X, Z = cbind(z1 = z1, z2 = z2))
}

# The skeleton a real fit builds, obtained cheaply, so the eligibility
# predicate can be asked about a state the driver would actually have made.
.cc_fit <- function(d, ...) suppressMessages(suppressWarnings(
  fit_mixture(indicators = d$X, n_classes = 2, measurement = "binary",
              predictors = d$Z, n_steps = 1, standard_errors = FALSE,
              order_by_size = FALSE, ...)))

test_that(".pattern_index() widens its key when given a second table", {
  X <- rbind(c(1, 0), c(1, 0), c(1, 0), c(0, 1))
  Z <- rbind(c(1),    c(1),    c(2),    c(1))

  # On the responses alone the first three rows are one pattern.
  expect_identical(.pattern_index(X)$idx, c(1L, 1L, 1L, 2L))
  # With the covariate in the key the third splits off, because its class
  # probabilities are not the same as the first two's.
  expect_identical(.pattern_index(X, Z)$idx, c(1L, 1L, 2L, 3L))
  # Z = NULL is the original key, byte for byte.
  expect_identical(.pattern_index(X, NULL), .pattern_index(X))
})

test_that(".collapse_patterns() collapses a discrete-covariate fit and declines the rest", {
  d  <- .cc_sim()
  st <- .cc_fit(d, n_init = 1, max_iter = 1L, refine = FALSE)

  coll <- .collapse_patterns(st, d$X, d$Z)
  expect_false(is.null(coll))
  expect_lt(nrow(coll$X), nrow(d$X) / 2)
  expect_equal(nrow(coll$Y), nrow(coll$X))
  expect_equal(sum(coll$w), nrow(d$X))
  # The pattern table rebuilds both tables exactly.
  pat <- .pattern_index(d$X, d$Z)
  expect_equal(coll$X[pat$idx, ], d$X)
  expect_equal(coll$Y[pat$idx, ], d$Z)

  # A missing covariate is refused: complete_covariates() imputes from an
  # unweighted mean, which the pattern table would silently change.
  Zna <- d$Z; Zna[1, 1] <- NA
  expect_null(.collapse_patterns(st, d$X, Zna))

  # Continuous covariates make almost every row unique, so the payoff gate
  # declines the collapse without needing a rule of its own.
  expect_null(.collapse_patterns(st, d$X, cbind(stats::rnorm(nrow(d$X)))))

  # No structural model of any other kind travels this path.
  other <- st; class(other$sm) <- "distal_continuous"
  expect_null(.collapse_patterns(other, d$X, d$Z))
})

test_that("the collapsed table gives the full sample's log-likelihood", {
  # The whole claim in one assertion: at the *same* parameters, the weighted
  # sum over the pattern table and the plain sum over the n cases are the same
  # number. Comparing two searches instead would compare their starting points.
  d  <- .cc_sim()
  st <- .cc_fit(d, n_init = 2, random_state = 4)

  coll <- .collapse_patterns(st, d$X, d$Z)
  expect_false(is.null(coll))

  full_ll <- sum(e_step(st, d$X, d$Z)$log_prob_norm)

  st_c <- st
  st_c$sample_weights <- coll$w
  coll_ll <- sum(coll$w * e_step(st_c, coll$X, coll$Y)$log_prob_norm)

  expect_equal(coll_ll, full_ll)
})

test_that("a collapsed covariate fit reports the log-likelihood it actually has", {
  # Self-consistency after .expand_patterns() puts the full sample back: the
  # number the fit prints must be the one its own parameters imply over all n
  # cases, not the one the pattern table produced.
  d   <- .cc_sim()
  fit <- .cc_fit(d, n_init = 4, random_state = 7)

  expect_equal(length(fit$sample_weights), nrow(d$X))
  expect_equal(fit$metrics$ll,
               sum(e_step(fit, d$X, d$Z)$log_prob_norm),
               tolerance = 1e-8)
})

# Ordinal items, and direct covariate effects on the items, take the same
# economy: the ordinal emission is the multinoulli one with a ragged category
# count, and a DIF emission's item probabilities are a function of the
# response pattern and the case's own item covariates, which join the key.
.cc_ordinal_sim <- function(seed = 5, n = 1500) {
  set.seed(seed)
  z  <- stats::rbinom(n, 1L, 0.5)
  s  <- 1L + stats::rbinom(n, 1L, stats::plogis(-0.4 + 0.9 * z))
  th <- rbind(c(-1.2, 1.0), c(0.8, 2.6))
  X  <- matrix(NA_integer_, n, 4)
  for (j in 1:4) {
    eta <- th[s, ] - if (j == 1) 0.8 * z else 0
    u   <- stats::runif(n)
    X[, j] <- 1L + (u < stats::plogis(eta[, 1])) + (u < stats::plogis(eta[, 2]))
  }
  colnames(X) <- paste0("i", 1:4)
  list(X = X, d = data.frame(z = z))
}

.cc_ordinal_fit <- function(d, dif = TRUE, ...) suppressMessages(suppressWarnings(
  fit_mixture(indicators = d$X, n_classes = 2, measurement = "ordinal",
              cats = rep(3L, 4), predictors = ~ z, data = d$d, n_steps = 1,
              predictors_items = if (dif) list(i1 = ~ z),
              standard_errors = FALSE, order_by_size = FALSE, ...)))

test_that("ordinal and DIF fits collapse, and the table gives the full log-likelihood", {
  d  <- .cc_ordinal_sim()
  Z  <- as.matrix(d$d)
  for (dif in c(FALSE, TRUE)) {
    st   <- .cc_ordinal_fit(d, dif = dif, n_init = 2, random_state = 3)
    coll <- .collapse_patterns(st, d$X, Z)
    expect_false(is.null(coll))
    expect_lt(nrow(coll$X), nrow(d$X) / 2)
    expect_equal(sum(coll$w), nrow(d$X))
    st_c <- st
    st_c$sample_weights <- coll$w
    if (dif) {
      # The covariate index is the representative rows' own.
      expect_equal(length(coll$dif$pat), nrow(coll$X))
      expect_equal(coll$dif$Z[, 1], coll$Y[, 1])
      expect_identical(coll$dif_orig, st$mm$dif)
      st_c$mm$dif <- coll$dif
    }
    expect_equal(sum(coll$w * e_step(st_c, coll$X, coll$Y)$log_prob_norm),
                 sum(e_step(st, d$X, Z)$log_prob_norm))
  }
})

test_that("a collapsed DIF fit hands back the full sample's covariate index", {
  d   <- .cc_ordinal_sim()
  fit <- .cc_ordinal_fit(d, n_init = 3, random_state = 9)
  expect_equal(length(fit$mm$dif$pat), nrow(d$X))
  expect_equal(length(fit$sample_weights), nrow(d$X))
  expect_equal(fit$metrics$ll,
               sum(e_step(fit, d$X, as.matrix(d$d))$log_prob_norm),
               tolerance = 1e-8)
})
