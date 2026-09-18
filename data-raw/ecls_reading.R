# Prepare the bundled `ecls_reading` dataset: five binary reading-proficiency
# indicators measured four times (fall and spring of kindergarten, fall and
# spring of first grade) on N = 3,575 children from the Early Childhood
# Longitudinal Study, Kindergarten Class of 1998-99 (ECLS-K), plus a poverty
# indicator.
#
# Source: U.S. Department of Education, National Center for Education
# Statistics (NCES), ECLS-K kindergarten-first grade public-use child file
# (https://nces.ed.gov/ecls/dataproducts.asp; no registration required).
# ECLS-K public-use data are U.S. federal government works and are in the
# public domain (17 U.S.C. sec. 105); the only obligation on users is not to
# attempt to identify respondents. The analytic extract in
# data-raw/ecls_reading/dp.analytic.dat defines the five indicators the way
# Kaplan (2008) does in his stage-sequential analysis of these data, and was
# obtained as supplementary material to a later published re-analysis of the
# same indicators; it is not regenerated here from the full public-use file.
#
# Layout of dp.analytic.dat: fixed-width, one column of width 1 (poverty) then
# twenty columns of width 2 (5 items x 4 occasions, occasion-major). A blank
# field is a fixed-format zero, NOT a missing value: there is no missing-value
# marker in this file, and reading blanks as zero reproduces the published
# univariate proportions exactly (every occasion's counts sum to 3,575).

raw_path <- "data-raw/ecls_reading/dp.analytic.dat"
d <- utils::read.fwf(raw_path, widths = c(1, rep(2, 20)))
d[is.na(d)] <- 0L

items <- c("letters", "beginning", "ending", "sight", "context")
ecls_reading <- data.frame(
  poverty = as.integer(d[, 1]),
  lapply(d[, -1], as.integer)
)
names(ecls_reading)[-1] <- paste0(rep(items, 4), "_t", rep(1:4, each = 5))

stopifnot(nrow(ecls_reading) == 3575L,
          all(ecls_reading$poverty %in% 0:1),
          all(as.matrix(ecls_reading[, -1]) %in% 0:1))

usethis::use_data(ecls_reading, overwrite = TRUE)
