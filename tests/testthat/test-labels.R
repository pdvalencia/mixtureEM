# Display-label shortening for printed tables (.shorten_labels and friends).

test_that(".shorten_labels leaves fitting labels untouched", {
  x <- c("age", "sexo.M", "a_name_of_28_characters_xxxx")
  out <- .shorten_labels(x, width = 28L)
  expect_identical(out, x)
  expect_null(attr(out, "legend"))
})

test_that(".shorten_labels abbreviates long labels and keeps a key", {
  x   <- c("age", "sexual_orientation.Not heterosexual")
  out <- .shorten_labels(x, width = 28L)

  expect_identical(out[1], "age")
  expect_lte(nchar(out[2]), 28L)
  legend <- attr(out, "legend")
  expect_identical(unname(legend), x[2])
  expect_identical(names(legend), out[2])
})

test_that(".shorten_labels disambiguates near-identical long names", {
  x <- c(paste0(strrep("very_long_prefix_", 3), "A"),
         paste0(strrep("very_long_prefix_", 3), "B"))
  out <- .shorten_labels(x, width = 20L)
  expect_false(out[1] == out[2])
  expect_identical(length(attr(out, "legend")), 2L)
})

test_that(".shorten_labels survives non-ASCII labels", {
  x   <- c(paste0("orientación_sexual.No heterosexual, en absoluto"))
  out <- expect_no_warning(.shorten_labels(x, width = 20L))
  expect_lte(nchar(out[1]), 20L)
})

test_that("summary tables shorten long dummy names but return full names", {
  set.seed(3)
  n  <- 200
  cl <- rbinom(n, 1, 0.4)
  items <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cl == 1, 0.85, 0.15)))
  colnames(items) <- paste0("item", 1:4)
  covs <- data.frame(
    sexual_orientation = factor(ifelse(
      rbinom(n, 1, plogis(cl - 0.4)) == 1, "Heterosexual",
      "Not heterosexual, including questioning")),
    age = rnorm(n)
  )

  fit  <- suppressMessages(fit_mixture(items, n_classes = 2,
                                       measurement = "binary",
                                       n_init = 3,
                                       random_state = 4))
  fitc <- suppressMessages(add_covariates(fit, covs))

  out_lines <- capture.output(sv <- summary(fitc))

  # The 51-character dummy name is not printed at full length anywhere,
  # and the key that decodes the abbreviation is.
  long_name <- "sexual_orientation.Not heterosexual, including questioning"
  expect_false(any(grepl(long_name, out_lines, fixed = TRUE) &
                     !grepl("Abbreviated names", out_lines, fixed = TRUE) &
                     !grepl(" = ", out_lines, fixed = TRUE)))
  expect_true(any(grepl("Abbreviated names:", out_lines, fixed = TRUE)))

  # The returned data frame carries the full, unshortened name.
  expect_true(long_name %in% sv$coefficients$term)
})

test_that("measurement_summary shortens very long indicator names", {
  set.seed(5)
  items <- matrix(rbinom(400, 1, 0.5), nrow = 100)
  colnames(items) <- paste0("an_extremely_verbose_indicator_name_number_", 1:4)
  fit <- suppressMessages(fit_mixture(items, n_classes = 2,
                                      measurement = "binary",
                                      n_init = 2,
                                      random_state = 6))
  out_lines <- capture.output(msdf <- measurement_summary(fit))
  expect_true(any(grepl("Abbreviated names:", out_lines, fixed = TRUE)))
  # Full names preserved in the returned data frame.
  expect_setequal(unique(msdf$item), colnames(items))
})
