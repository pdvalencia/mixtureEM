# print.mixture_comparison(): the short default view, the full view, and the
# formatting of BF and cmP. Built from a hand-made table so no model is fitted.

.cp_obj <- function(unrep = c(FALSE, FALSE, FALSE), vlmr = FALSE) {
  tab <- data.frame(Classes = 1:3, LL = c(-500, -450, -440), Params = c(5, 11, 17),
                    AIC = c(1010, 922, 914), BIC = c(1030, 960, 980),
                    CAIC = c(1035, 971, 997), AWE = c(1070, 1040, 1100),
                    AIC3 = c(1015, 933, 931), ICL = c(1030, 990, 1020),
                    SABIC = c(1014, 925, 926), Entropy = c(1, 0.81, 0.77),
                    Unreplicated = unrep)
  tab <- .bf_cmp(tab)
  if (vlmr) { tab$VLMR_LR <- c(100, 20, NA); tab$VLMR_p <- c(0.0001, 0.2, NA) }
  out <- list(fit_table = tab, models = list(), best_k = 2L)
  class(out) <- "mixture_comparison"
  out
}

test_that("the default view keeps the short set of columns", {
  out <- capture.output(r <- withVisible(print(.cp_obj())))
  hdr <- out[grep("Classes", out)]
  for (nm in c("LL", "Params", "BIC", "SABIC", "AIC3", "Entropy"))
    expect_match(hdr, paste0("\\b", nm, "\\b"))
  for (nm in c("AIC ", "CAIC", "AWE", "ICL", "BF", "cmP"))
    expect_false(grepl(nm, paste0(hdr, " "), fixed = TRUE))
  expect_true(any(grepl("Best model according to BIC: 2 classes", out, fixed = TRUE)))
  expect_true(any(grepl("columns = \"all\"", out, fixed = TRUE)))
  expect_false(r$visible)
})

test_that("the full view adds every criterion and caps BF at its reading bands", {
  x   <- .cp_obj()
  # BF row 1 = exp(-(1030 - 960) / 2) ~ 6e-16; row 2 = exp((980 - 960) / 2) ~ 22026.
  out <- capture.output(print(x, columns = "all"))
  txt <- paste(out, collapse = "\n")
  for (nm in c("AIC", "CAIC", "AWE", "ICL", "BF", "cmP"))
    expect_match(txt, paste0("\\b", nm, "\\b"))
  expect_match(txt, "< 0.1", fixed = TRUE)
  expect_match(txt, "> 100", fixed = TRUE)
  expect_match(txt, "< 0.001", fixed = TRUE)
  expect_false(grepl("e+", txt, fixed = TRUE))
  expect_false(grepl("e-", txt, fixed = TRUE))
  # Between the bands BF prints to two decimals.
  x$fit_table$BF[1] <- 3.456
  expect_match(paste(capture.output(print(x, columns = "all")), collapse = "\n"),
               "3.46", fixed = TRUE)
})

test_that("the table itself keeps full precision", {
  x <- .cp_obj()
  before <- x$fit_table
  capture.output(print(x, columns = "all"))
  expect_identical(x$fit_table, before)
  expect_lt(x$fit_table$BF[1], 1e-10)
})

test_that("VLMR columns are shown whenever they were computed", {
  out <- paste(capture.output(print(.cp_obj(vlmr = TRUE))), collapse = "\n")
  expect_match(out, "VLMR_p", fixed = TRUE)
  expect_match(out, "0.000", fixed = TRUE)
})

test_that("a longitudinal comparison names its unit", {
  x <- .cp_obj()
  attr(x, "unit") <- "statuses"
  expect_true(any(grepl("BIC: 2 statuses", capture.output(print(x)), fixed = TRUE)))
})
