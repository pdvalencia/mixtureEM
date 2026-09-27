# Every model the package can fit is either taken to the maximum by the
# Newton-type finish (R/qn_finish.R, R/qn_pack.R) or named in .qn_refused with
# the reason it is not. This file is the guard on that: it reads the
# descriptors build_emission() accepts straight from its source, so a new
# measurement or structural model fails here until it is given a packing or a
# written refusal.
#
# For each model it checks the three things the finish relies on: that the fit
# is packable, that packing and unpacking give back the same parameters, and
# that the packed model reproduces the fit's own log-likelihood.

.cov_q <- function(e) suppressMessages(suppressWarnings(e))

.cov_bin <- function(n = 240, J = 4, na = FALSE, seed = 1) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE)
  P  <- rbind(rep(.8, J), rep(.2, J))
  X  <- t(vapply(cl, function(k) rbinom(J, 1, P[k, ]), numeric(J)))
  if (na) X[sample(length(X), 20)] <- NA
  colnames(X) <- paste0("u", seq_len(J))
  X
}
.cov_cont <- function(n = 240, J = 3, na = FALSE, seed = 2) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE)
  X  <- matrix(rnorm(n * J), n) + 2 * (cl == 2)
  if (na) X[sample(length(X), 20)] <- NA
  X
}
.cov_count <- function(n = 240, J = 3, na = FALSE, seed = 3) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE)
  X  <- matrix(rpois(n * J, ifelse(cl == 1, 0.5, 3)), n)
  if (na) X[sample(length(X), 20)] <- NA
  X
}
.cov_cat <- function(n = 240, J = 3, na = FALSE, seed = 4) {
  set.seed(seed)
  cl <- sample(1:2, n, TRUE)
  X  <- t(vapply(cl, function(k) vapply(seq_len(J), function(j)
    sample.int(3, 1, prob = if (k == 1) c(.7, .2, .1) else c(.1, .3, .6)),
    integer(1)), integer(J)))
  if (na) X[sample(length(X), 20)] <- NA
  X
}
.cov_fm <- function(X, ...) .cov_q(fit_mixture(X, n_classes = 2, n_init = 2,
                                               random_state = 1, n_cores = 1, ...))

# The objective at the fit's own packed point against the fit's own
# log-likelihood plus prior, and the round trip.
.cov_check_mixture <- function(fit, Y, label) {
  expect_true(.qn_mixture_supported(fit, Y), label = label)
  pr <- .qn_mixture_problem(fit, fit$data, Y)
  st <- pr$unpack(pr$par)
  expect_equal(.qn_mixture_layout(st, Y)$par, pr$par, tolerance = 1e-8,
               label = paste(label, "round trip"))
  ll <- sum(fit$sample_weights * e_step(fit, fit$data, Y)$log_prob_norm)
  expect_equal(sum(pr$w * pr$case_ll(pr$par)), ll, tolerance = 1e-8,
               label = paste(label, "log-likelihood"))
}

# One small fit per class the finish must reach, and the structural data the
# finish sees for it (NULL on a stepwise fit, whose finish is step 1's).
.cov_fixtures <- list(
  bernoulli         = function() list(fit = .cov_fm(.cov_bin(), measurement = "binary")),
  bernoulli_nan     = function() list(fit = .cov_fm(.cov_bin(na = TRUE), measurement = "binary")),
  multinoulli       = function() list(fit = .cov_fm(.cov_cat(), measurement = "categorical")),
  multinoulli_nan   = function() list(fit = .cov_fm(.cov_cat(na = TRUE), measurement = "categorical")),
  gaussian_unit     = function() list(fit = .cov_fm(.cov_cont(), measurement = "gaussian")),
  gaussian_unit_nan = function() list(fit = .cov_fm(.cov_cont(na = TRUE), measurement = "gaussian")),
  gaussian_diag     = function() list(fit = .cov_fm(.cov_cont(), measurement = "continuous")),
  gaussian_diag_nan = function() list(fit = .cov_fm(.cov_cont(na = TRUE), measurement = "continuous")),
  poisson           = function() list(fit = .cov_fm(.cov_count(), measurement = "count")),
  poisson_nan       = function() list(fit = .cov_fm(.cov_count(na = TRUE), measurement = "count")),
  nested            = function() list(fit = .cov_fm(cbind(.cov_bin(J = 3), .cov_cont(J = 2)),
                                        measurement = list(binary = 1:3, continuous = 4:5))),
  group_blocks      = function() {
    X <- .cov_cont(); g <- factor(rep(1:2, length.out = nrow(X)))
    list(fit = .cov_fm(X, measurement = "continuous", group = g,
                       group_effects = "measurement", group_invariant_params = "means"))
  },
  lcga              = function() {
    set.seed(5); Y <- matrix(rbinom(200 * 4, 1, rep(c(.2, .7), each = 100)), 200)
    list(fit = .cov_q(fit_lcga(Y, times = 4, n_classes = 2, family = "binomial",
                               n_init = 2, random_state = 1)))
  },
  structured_normal = function() {
    set.seed(6); Y <- outer(rnorm(200, rep(c(0, 3), each = 100)), 0:3 * .5, "+") +
      matrix(rnorm(800, sd = .5), 200)
    list(fit = .cov_q(fit_gmm(Y, times = 4, n_classes = 2, n_init = 2,
                              random_state = 1)))
  },
  bernoulli_dif     = function() {
    X <- .cov_bin(); d <- data.frame(z = rbinom(nrow(X), 1, .5))
    list(fit = .cov_fm(X, measurement = "binary", predictors = ~ z, data = d,
                       predictors_items = list(u1 = ~ z)), one_step = TRUE)
  },
  covariate         = function() {
    X <- .cov_bin()
    list(fit = .cov_fm(X, measurement = "binary", predictors = rnorm(nrow(X)),
                       n_steps = 1), one_step = TRUE)
  },
  group_prevalence  = function() {
    X <- .cov_bin(); g <- factor(rep(1:2, length.out = nrow(X)))
    list(fit = .cov_fm(X, measurement = "binary", group = g,
                       group_effects = "prevalence", group_prevalence_equal = 1,
                       n_steps = 1), one_step = TRUE)
  },
  distal_continuous = function() {
    X <- .cov_bin()
    list(fit = .cov_fm(X, measurement = "binary", outcome = rnorm(nrow(X)),
                       n_steps = 1), one_step = TRUE)
  },
  distal_continuous_pooled = function() {
    X <- .cov_bin(); z <- data.frame(z = rnorm(nrow(X)))
    list(fit = .cov_fm(X, measurement = "binary", outcome = rnorm(nrow(X)),
                       outcome_covariates = z, n_steps = 1), one_step = TRUE)
  },
  # `slopes = "class_specific"` now fits distal_continuous_pooled, so this
  # engine is reached through its descriptor only.
  distal_continuous_regression = function() {
    X <- .cov_bin()
    Y <- cbind(y = rnorm(nrow(X)), z = rnorm(nrow(X)))
    list(fit = .cov_q(fit_mixture_internal(
           X, Y, n_components = 2, measurement = "binary",
           structural = "distal_continuous_regression", n_steps = 1,
           n_init = 2, random_state = 1, n_cores = 1)), one_step = TRUE)
  },
  distal_categorical = function() {
    X <- .cov_bin()
    list(fit = .cov_fm(X, measurement = "binary",
                       outcome = factor(sample(1:3, nrow(X), TRUE)), n_steps = 1),
         one_step = TRUE)
  },
  distal_pooled     = function() {
    X <- .cov_bin(); z <- data.frame(z = rnorm(nrow(X)))
    list(fit = .cov_fm(X, measurement = "binary",
                       outcome = factor(sample(1:2, nrow(X), TRUE)),
                       outcome_covariates = z, n_steps = 1), one_step = TRUE)
  },
  distal_regression = function() {
    X <- .cov_bin(); z <- data.frame(z = rnorm(nrow(X)))
    list(fit = .cov_fm(X, measurement = "binary",
                       outcome = factor(sample(1:2, nrow(X), TRUE)),
                       outcome_covariates = z, slopes = "class_specific",
                       n_steps = 1), one_step = TRUE)
  }
)

# Descriptor -> the class it builds. Every string build_emission() tests
# `descriptor` against must be here. Classes reached only through fit_lta()
# (ordinal, time_blocks) are checked in the LTA test below.
.cov_classes <- c(
  bernoulli = "bernoulli", binary = "bernoulli",
  bernoulli_nan = "bernoulli_nan", binary_nan = "bernoulli_nan",
  multinoulli = "multinoulli", categorical = "multinoulli",
  multinoulli_nan = "multinoulli_nan", categorical_nan = "multinoulli_nan",
  gaussian_unit = "gaussian_unit", gaussian = "gaussian_unit",
  gaussian_unit_nan = "gaussian_unit_nan", gaussian_nan = "gaussian_unit_nan",
  gaussian_diag = "gaussian_diag", continuous = "gaussian_diag",
  gaussian_diag_nan = "gaussian_diag_nan", continuous_nan = "gaussian_diag_nan",
  poisson = "poisson", count = "poisson",
  poisson_nan = "poisson_nan", count_nan = "poisson_nan",
  ordinal = "lta:ordinal", ordinal_nan = "lta:ordinal",
  time_blocks = "lta:time_blocks", group_blocks = "group_blocks",
  lcga = "lcga", structured_normal = "structured_normal",
  gmm = "structured_normal",
  covariate = "covariate", predict_class = "covariate",
  distal_regression = "distal_regression",
  categorical_outcome_moderated = "distal_regression",
  distal_pooled = "distal_pooled", categorical_outcome_adjusted = "distal_pooled",
  categorical_outcome = "distal_categorical",
  distal_continuous = "distal_continuous", continuous_outcome = "distal_continuous",
  distal_continuous_regression = "distal_continuous_regression",
  continuous_outcome_moderated = "distal_continuous_regression",
  distal_continuous_pooled = "distal_continuous_pooled",
  continuous_outcome_adjusted = "distal_continuous_pooled",
  group_prevalence = "group_prevalence")

test_that("every descriptor build_emission() accepts is accounted for", {
  txt  <- paste(deparse(body(build_emission)), collapse = "\n")
  lits <- regmatches(txt, gregexpr(
    'descriptor (%in%|==) (c\\([^)]*\\)|"[^"]*")', txt))[[1]]
  desc <- unique(gsub('"', "", unlist(regmatches(lits, gregexpr('"[^"]*"', lits)))))
  expect_setequal(desc, names(.cov_classes))
  cls <- setdiff(unique(.cov_classes), c("lta:ordinal", "lta:time_blocks"))
  # A class is covered by a fixture below, or refused with a reason.
  covered <- c(names(.cov_fixtures), names(.qn_refused))
  expect_true(all(cls %in% covered),
              label = paste("uncovered:", paste(setdiff(cls, covered), collapse = ", ")))
  expect_true(all(nzchar(.qn_refused)))
})

test_that("every mixture-engine model is packable and packs losslessly", {
  for (nm in names(.cov_fixtures)) {
    fx  <- .cov_fixtures[[nm]]()
    fit <- fx$fit
    Y   <- if (isTRUE(fx$one_step)) fit$Y else NULL
    .cov_check_mixture(fit, Y, nm)
    # The fixture exercises the class it is named for.
    classes <- c(class(fit$mm), class(fit$sm),
                 unlist(lapply(fit$mm$models, class)))
    expect_true(nm %in% classes, label = paste(nm, "class present"))
  }
})

.cov_lta_data <- function(n = 200, Tn = 2, J = 3, gen) {
  set.seed(9)
  s <- matrix(sample(1:2, n * Tn, TRUE), n)
  X <- NULL
  for (t in seq_len(Tn)) X <- cbind(X, vapply(seq_len(J), function(j) gen(s[, t]),
                                               numeric(n)))
  X
}
.cov_lta <- function(X, times = 2, ...)
  .cov_q(fit_lta(X, n_statuses = 2, times = times, n_init = 2, random_state = 1,
                 n_cores = 1, standard_errors = FALSE, ...))

test_that("every latent transition model is finished, or refused by name", {
  bin  <- .cov_lta_data(gen = function(s) rbinom(length(s), 1, c(.2, .8)[s]))
  # A mixture over chains is identified from three occasions on.
  bin3 <- .cov_lta_data(Tn = 3, gen = function(s) rbinom(length(s), 1, c(.2, .8)[s]))
  cont <- .cov_lta_data(gen = function(s) rnorm(length(s), 2 * s))
  cnt  <- .cov_lta_data(gen = function(s) rpois(length(s), c(.5, 3)[s]))
  cat3 <- .cov_lta_data(gen = function(s) vapply(s, function(k)
    sample.int(3, 1, prob = if (k == 1) c(.7, .2, .1) else c(.1, .2, .7)), integer(1)))
  fits <- list(
    binary      = .cov_lta(bin, measurement = "binary"),
    ordinal     = .cov_lta(cat3, measurement = "ordinal"),
    continuous  = .cov_lta(cont, measurement = "continuous"),
    categorical = .cov_lta(cat3, measurement = "categorical"),
    count       = .cov_lta(cnt, measurement = "count"),
    two_chains  = .cov_lta(bin3, times = 3, measurement = "binary", n_classes = 2),
    mover_stayer = .cov_lta(bin3, times = 3, measurement = "binary",
                            mover_stayer = TRUE))
  for (nm in names(fits))
    expect_true(.qn_lta_supported(fits[[nm]]), label = nm)
  # The packings the finite-difference route uses, round-tripped.
  for (nm in c("ordinal", "categorical", "count")) {
    mm <- fits[[nm]]$mm
    pk <- .qn_pack_mm(mm)
    expect_false(is.null(pk), label = nm)
    back <- .qn_unpack_mm(mm, pk$par)
    expect_equal(.qn_pack_mm(back)$par, pk$par, tolerance = 1e-10, label = nm)
    expect_equal(.lta_emission_loglik(back, fits[[nm]]$data),
                 .lta_emission_loglik(mm, fits[[nm]]$data), tolerance = 1e-10,
                 label = nm)
  }
  # A freely loading random intercept is finished on binary items and is the
  # refusal on the list on ordinal ones.
  for (kind in c("binary", "continuous")) {
    ri <- .cov_lta(bin, measurement = "binary", random_intercept = kind)
    expect_true(.qn_lta_supported(ri), label = paste("binary items,", kind, "RI"))
  }
  ri_ord <- .cov_lta(cat3, measurement = "ordinal", random_intercept = "continuous")
  expect_false(.qn_lta_supported(ri_ord))
  expect_true("lta_random_intercept" %in% names(.qn_refused))
})
