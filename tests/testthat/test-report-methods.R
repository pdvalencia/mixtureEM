# report_methods() and grolts_report() print what the fit already stores. The
# snapshots guard the wording with every digit masked, so they do not move with
# a platform's last decimal; the numbers are guarded by the traceability test,
# which requires every number in the text to be one the object holds.

mask_digits <- function(x) gsub("[0-9]+", "#", x)

lca_fixture <- function() {
  set.seed(11)
  z <- rbinom(150, 1, 0.4)
  p <- ifelse(z == 1, 0.85, 0.15)
  X <- matrix(rbinom(150 * 6, 1, rep(p, 6)), ncol = 6)
  X[1, ] <- NA
  X[2:6, 2] <- NA
  suppressWarnings(fit_mixture(X, n_classes = 2, measurement = "binary",
                               n_init = 12, random_state = 4))
}

test_that("every number in the LCA report is one the fit stores", {
  fit <- lca_fixture()
  rep <- report_methods(fit)
  expect_s3_class(rep, "methods_report")

  m   <- fit$metrics
  cs  <- class_sizes(fit)
  utils::capture.output(cd <- classification_diagnostics(fit))
  pp  <- diag(cd$ave_pp)
  md  <- fit$missing_data
  ver <- as.character(utils::packageVersion("mixtureEM"))
  stored <- c(
    sprintf("%.2f", c(m$ll, m$aic, m$bic, m$sabic)),
    sprintf("%.3f", c(m$entropy, min(pp), max(pp))),
    sprintf("%.1f%%", 100 * min(fit$weights)),
    sprintf("%.0f", min(cs$n_expected)),
    sprintf("%.1f%%", 100 * md$prop_missing),
    as.character(c(m$n_params, m$n_starts, m$n_replicated, fit$n_iter,
                   fit$n_components, md$n_input_rows, md$n_empty_rows,
                   md$n_input_rows - md$n_empty_rows, md$n_missing,
                   md$n_cells, md$n_items_affected, ncol(fit$data),
                   fit$random_state, unlist(fit$bayes_constants))),
    ver, paste(R.version$major, R.version$minor, sep = "."),
    "1977", "2013", "2006", "2008")
  text <- paste(rep[-length(rep)], collapse = " ")  # the call line holds arguments
  nums <- regmatches(text, gregexpr("-?[0-9]+(\\.[0-9]+)*%?", text))[[1]]
  # Author-year citations and the package version split into parts are fine.
  stray <- setdiff(nums, c(stored, unlist(strsplit(ver, ".", fixed = TRUE))))
  expect_identical(stray, character(0))
})

test_that("the report states the starts, the prior and the missing data", {
  fit  <- lca_fixture()
  text <- paste(report_methods(fit), collapse = " ")
  expect_match(text, "of 150 supplied; 1 case had no observed indicator")
  expect_match(text, "The search used 12 random starts")
  expect_match(text, "full-information maximum likelihood")
  expect_match(text, "1 pseudo-observation for the class proportions")
  expect_match(text, "random_state = 4")
  # Never graded: no thresholds, no verdicts.
  expect_no_match(text, "adequate|acceptable|good|poor|sufficient", ignore.case = TRUE)
  expect_true(any(grepl("^Masyn", attr(report_methods(fit), "references"))))
})

test_that("a plain-ML fit says so and an unseeded fit says it cannot be reproduced", {
  set.seed(3)
  X <- matrix(rbinom(400, 1, 0.5), ncol = 4)
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2, measurement = "binary",
    n_init = 3, bayes_constants = list(latent = 0, categorical = 0)))
  text <- paste(report_methods(fit), collapse = " ")
  expect_match(text, "with no prior on any parameter")
  expect_match(text, "No random_state was set")
})

test_that("the report wording is stable for each model family", {
  fit <- lca_fixture()
  expect_snapshot(print(report_methods(fit)), transform = mask_digits)

  set.seed(5)
  y <- matrix(rnorm(400, sd = 0.7), 100) +
    outer(rep(c(0, 2), each = 50), 0:3, function(a, t) a + (a / 2) * t)
  lcga <- suppressWarnings(fit_lcga(y, times = 4, n_classes = 2,
    family = "gaussian", n_init = 3, random_state = 1))
  expect_snapshot(print(report_methods(lcga)), transform = mask_digits)
  expect_snapshot(print(grolts_report(lcga)), transform = mask_digits)

  set.seed(22)
  Xl <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  lta <- suppressWarnings(fit_lta(Xl, n_statuses = 2, times = 3,
    measurement = "binary", n_init = 2, random_state = 2,
    standard_errors = FALSE))
  expect_snapshot(print(report_methods(lta)), transform = mask_digits)
})

test_that("a comparison adds the enumeration paragraph and the call that rebuilds it", {
  set.seed(8)
  X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
  utils::capture.output(cmp <- suppressWarnings(compare_mixtures(X, k_range = 1:2,
    measurement = "binary", n_init = 3, random_state = 9)))
  text <- paste(report_methods(cmp), collapse = " ")
  expect_match(text, "Models with 1 to 2 classes were compared, including the one-class model")
  expect_match(text, sprintf("BIC was lowest for K = %d, the model reported here", cmp$best_k))
  expect_match(text, "The models are refitted by: compare_mixtures(", fixed = TRUE)
  expect_match(text, "random_state = 9")
  expect_error(report_methods(cmp, k = 5), "holds no K = 5 model")
})

test_that("a three-step fit names its correction by its literature name", {
  set.seed(2)
  z <- rbinom(200, 1, 0.5)
  X <- matrix(rbinom(200 * 5, 1, rep(ifelse(z == 1, 0.8, 0.2), 5)), ncol = 5)
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2, measurement = "binary",
    predictors = matrix(rnorm(200) + z, ncol = 1), n_steps = 3,
    correction = "ML", n_init = 3, random_state = 1))
  rep <- report_methods(fit)
  text <- paste(rep, collapse = " ")
  expect_match(text, "maximum-likelihood correction (Vermunt, 2010)", fixed = TRUE)
  expect_match(text, "class membership on the predictors")
  # The fit indices are the measurement model's, never step 3's.
  expect_match(text, sprintf("%.2f", fit$step1_metrics$ll), fixed = TRUE)
  # How many cases step 3 moved, read off classification_diagnostics().
  capture.output(cd <- classification_diagnostics(fit))
  expect_match(text, sprintf("The modal class of %d case", round(cd$class_shift$n_moved)),
               fixed = TRUE)
  expect_no_match(text, "Mahalanobis")   # binary indicators
})

test_that("a continuous fit reports its smallest class separation", {
  set.seed(31)
  z <- sample(1:2, 200, TRUE)
  X <- outer(c(0, 2)[z], c(1, 1, 1)) + matrix(rnorm(600), 200)
  fit <- suppressWarnings(fit_mixture(X, n_classes = 2,
    measurement = "continuous", n_init = 3, random_state = 1))
  capture.output(cd <- classification_diagnostics(fit))
  text <- paste(report_methods(fit), collapse = " ")
  expect_match(text, sprintf("two class means, against the pooled within-class variance, was %.2f",
                             attr(cd$separation, "min")), fixed = TRUE)
})

test_that("grolts_report() walks all 21 items, prints no score, and refuses a non-trajectory fit", {
  set.seed(5)
  y <- matrix(rnorm(400, sd = 0.7), 100) +
    outer(rep(c(0, 2), each = 50), 0:3, function(a, t) a + (a / 2) * t)
  fit <- suppressWarnings(fit_lcga(y, times = 4, n_classes = 2,
    family = "gaussian", n_init = 3, random_state = 1))
  g <- grolts_report(fit)
  expect_s3_class(g, "grolts_report")
  expect_identical(nrow(g), 21L)
  expect_setequal(unique(g$status), c("from the fit", "you supply", "not applicable"))
  out <- paste(utils::capture.output(print(g)), collapse = " ")
  expect_no_match(out, "[0-9]+ (of|out of) 21")
  expect_identical(g$status[g$item == "8"], "not applicable")
  expect_match(g$answer[g$item == "1"], "Time scores 0, 1, 2, 3")

  expect_error(grolts_report(lca_fixture()), "for trajectory models")
})
