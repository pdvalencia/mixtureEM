# Reproduces Ventura-Leon et al. (2025), Tables 2-3: a 4-class LCA of the
# bundled `ventura_leon` data should recover their published class sizes,
# item-response probabilities, and parameter count. This is both a check on
# the data-raw/ventura_leon.R recode and a regression test for the core EM
# engine against a real, independently published solution.

test_that("4-class LCA of ventura_leon matches Ventura-Leon et al. (2025) Tables 2-3", {
  skip_on_cran()  # 400 cases x 30 random starts is slow for a CRAN check

  # The 16 infidelity items, in the instrument's item1-item16 order (the
  # bundled data carries informative column names; positions 7-22).
  items <- ventura_leon[, 7:22]

  set.seed(1)
  fit <- fit_mixture(items, n_classes = 4, measurement = "binary",
                      n_init = 30, max_iter = 2000)

  expect_true(fit$converged)
  # Table 2: 4-class model has 67 free parameters (16 items x 4 classes + 3
  # class-weight parameters).
  expect_equal(fit$metrics$n_params, 67L)

  # Table 3 class prevalences, largest to smallest: Fidelity .425, Affective
  # interest .272, Infidelity .157, Sexual desire .144. Our EM consistently
  # converges (checked across several seeds) to a global optimum whose
  # Sexual-desire/Fidelity split differs from the published one by ~1pp,
  # most likely reflecting a distinct near-tied optimum in the original
  # multilevLCA fit rather than an error in either implementation.
  paper_prevalence <- c(.425, .272, .157, .144)
  expect_equal(sort(fit$weights, decreasing = TRUE), paper_prevalence,
               tolerance = 0.03)

  # Table 3 item-response ("Yes"/endorsement) probabilities, in the same
  # size-sorted class order (Fidelity, Affective interest, Infidelity,
  # Sexual desire).
  paper_pis <- rbind(
    item1  = c(.15, .62, .98, .76),
    item2  = c(.18, .42, .94, .25),
    item3  = c(.17, .56, .92, .51),
    item4  = c(.09, .43, .92, .39),
    item5  = c(.11, .39, .78, .20),
    item6  = c(.06, .43, .81, .36),
    item7  = c(.23, .75, .97, .85),
    item8  = c(.06, .68, .95, .83),
    item9  = c(.00, .13, 1.00, .02),
    item10 = c(.00, .13, 1.00, .06),
    item11 = c(.00, .03, .90, .85),
    item12 = c(.01, .00, .92, .85),
    item13 = c(.00, .03, .90, .65),
    item14 = c(.15, .67, 1.00, .98),
    item15 = c(.00, .08, .87, .00),
    item16 = c(.01, .07, .89, .83)
  )

  ord <- order(fit$weights, decreasing = TRUE)
  pis <- fit$mm$parameters$pis[ord, ]
  colnames(pis) <- colnames(items)
  expect_equal(unname(t(pis)), unname(paper_pis), tolerance = 0.06)
})
