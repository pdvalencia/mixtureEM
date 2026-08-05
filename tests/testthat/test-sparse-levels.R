# Factor handling in prepare_covariates(): unused levels are dropped, sparse
# levels are flagged.

test_that("unused factor levels produce no dummy column", {
  gender <- factor(rep(c("Male", "Female"), 50),
                   levels = c("Male", "Female", "Other"))
  # No "Other" cases at all — the classic subset-without-droplevels situation.
  expect_message(
    Y <- prepare_covariates(data.frame(gender = gender)),
    "dropped 1 unused level")
  expect_identical(colnames(Y), "gender.Female")
  expect_identical(attr(Y, "covariate_terms"), "gender")
})

test_that("sparse observed levels keep their dummy but warn", {
  gender <- factor(c(rep("Male", 60), rep("Female", 58), "Other", "Other"))
  expect_warning(
    Y <- prepare_covariates(data.frame(gender = gender)),
    "'Other' \\(n = 2\\)")
  expect_true("gender.Other" %in% colnames(Y))
})

test_that("well-populated factors pass through silently", {
  gender <- factor(rep(c("Male", "Female", "Other"), each = 20))
  expect_no_warning(expect_no_message(
    Y <- prepare_covariates(data.frame(gender = gender))))
  expect_identical(ncol(Y), 2L)
})
