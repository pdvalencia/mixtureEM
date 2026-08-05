# Reproduces Collins & Lanza (2010), sec. 2.4 and Table 4.5/2.7: a 5-class
# LCA of the bundled `yrbs2005` data should recover their published class
# sizes, item-response probabilities, and parameter count. This is both a
# check on the data-raw/yrbs2005.R recode and a regression test for the core
# EM engine against a real, independently published solution.

test_that("5-class LCA of yrbs2005 matches Collins & Lanza (2010) Table 2.7/4.5", {
  skip_on_cran()  # ~13,900 cases x 20 random starts is slow for a CRAN check

  items <- yrbs2005[, c("smoked_before_13", "smoked_daily_30d", "drove_drinking",
                        "first_drink_before_13", "binge_drink_30d",
                        "marijuana_before_13", "cocaine_ever", "glue_ever",
                        "meth_ever", "ecstasy_ever", "sex_before_13",
                        "sex_4plus_partners")]

  set.seed(1)
  fit <- fit_mixture(items, n_classes = 5, measurement = "binary",
                      n_init = 20, max_iter = 2000)

  expect_true(fit$converged)
  # Table 4.5: 5-class model, 54 parameters for the measurement model.
  expect_equal(fit$metrics$n_params, 64L)

  # Table 2.7 class prevalences, largest to smallest: Low Risk .67, Binge
  # Drinkers .14, Early Experimenters .09, High Risk .05, Sexual Risk-Takers
  # .04. `order_by_size = TRUE` (the default) sorts our classes the same way.
  book_prevalence <- c(.67, .14, .09, .05, .04)
  expect_equal(sort(fit$weights, decreasing = TRUE), book_prevalence,
               tolerance = 0.02)

  # Table 2.7 item-response probabilities in the same size-sorted class
  # order (Low Risk, Binge Drinkers, Early Experimenters, High Risk, Sexual
  # Risk-Takers), probability of a "Yes" response.
  book_pis <- rbind(
    smoked_before_13      = c(.04, .11, .76, .64, .17),
    smoked_daily_30d      = c(.02, .27, .31, .66, .12),
    drove_drinking        = c(.01, .42, .15, .45, .11),
    first_drink_before_13 = c(.14, .21, .79, .68, .39),
    binge_drink_30d       = c(.08, .74, .48, .79, .16),
    marijuana_before_13   = c(.01, .03, .46, .55, .22),
    cocaine_ever          = c(.00, .19, .07, .88, .03),
    glue_ever             = c(.06, .19, .22, .58, .04),
    meth_ever             = c(.00, .10, .02, .73, .01),
    ecstasy_ever          = c(.00, .11, .06, .64, .06),
    sex_before_13         = c(.01, .00, .18, .30, .81),
    sex_4plus_partners    = c(.06, .29, .24, .56, .83)
  )

  pis <- fit$mm$parameters$pis[order(fit$weights, decreasing = TRUE), ]
  colnames(pis) <- colnames(items)
  expect_equal(unname(t(pis)), unname(book_pis), tolerance = 0.02)
})
