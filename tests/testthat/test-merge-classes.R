# merge_classes(): post-processing on the posteriors, so the fit is untouched
# and the merged posteriors are sums of the fitted columns.

local({
  set.seed(11)
  X <- rbind(matrix(rnorm(300, 0), ncol = 3), matrix(rnorm(300, 3), ncol = 3))
  fit <- fit_mixture(X, n_classes = 3, measurement = "continuous",
                     n_init = 2, random_state = 1)

  test_that("merging leaves the fit untouched", {
    before <- fit
    m <- merge_classes(fit)
    expect_identical(fit, before)
    expect_s3_class(m, "class_merge")
    expect_identical(m$path$K, 3:1)
  })

  test_that("merged posteriors are sums of fitted columns and end at one", {
    m <- merge_classes(fit)
    resp <- exp(fit$log_resp)
    for (r in seq_along(m$posterior)) {
      expect_equal(rowSums(m$posterior[[r]]), rep(1, nrow(X)),
                   tolerance = 1e-12)
      for (j in seq_along(m$clusters[[r]]))
        expect_equal(m$posterior[[r]][, j],
                     rowSums(resp[, m$clusters[[r]][[j]], drop = FALSE]))
    }
    expect_equal(unname(m$posterior[[3]][, 1]), rep(1, nrow(X)))
    expect_equal(m$path$merge_entropy[3], 0)
  })

  test_that("the entropy at the fitted K is (ICL - BIC) / 2", {
    m <- merge_classes(fit)
    expect_equal(m$path$merge_entropy[1],
                 (fit$metrics$icl - fit$metrics$bic) / 2, tolerance = 1e-8)
    expect_true(all(diff(m$path$merge_entropy) <= 1e-12))
  })

  test_that("print and plot run", {
    m <- merge_classes(fit)
    expect_output(print(m), "Merge path")
    pdf(NULL); on.exit(dev.off())
    expect_invisible(plot(m))
  })
})

test_that("a one-class fit and an LTA fit are refused", {
  X <- matrix(rnorm(100), ncol = 2)
  f1 <- fit_mixture(X, n_classes = 1, measurement = "continuous")
  expect_error(merge_classes(f1), "nothing to merge")
  expect_error(merge_classes(structure(list(), class = "lta_model")),
               "latent statuses")
})
