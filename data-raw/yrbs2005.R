# Prepare the bundled `yrbs2005` dataset from the CDC's 2005 national Youth
# Risk Behavior Survey public-use file.
#
# Source: Centers for Disease Control and Prevention (2006). Youth Risk
# Behavior Survey, 2005 [national public-use data set]. Downloaded from
# https://ftp.cdc.gov/pub/data/yrbs/2005/ (no registration required; YRBS
# national data are U.S. federal government works and are in the public
# domain, 17 U.S.C. sec. 105).
#
# Column positions below are taken verbatim from
# data-raw/yrbs2005/YRBS_2005_SAS_Input_Program.sas (the CDC-supplied SAS
# input program shipped alongside the ASCII file). QN* variables are the
# CDC's own precomputed dichotomous recodes of the raw questionnaire items;
# using them (rather than recoding Q* ourselves) matches CDC's own analytic
# variables, which is also what Collins & Lanza (2010), sec. 1.8.1.2, 2.4,
# 4.6.4 and 5.12.2, used for their 12-item health-risk-behavior illustration.

raw_path <- "data-raw/yrbs2005/yrbs2005.dat"
lines <- readLines(raw_path)

get_col <- function(lines, start, width = 1L) {
  trimws(substr(lines, start, start + width - 1L))
}

to_binary <- function(x) {
  # CDC coding: "1" = yes, "2" = no, "." (or blank) = missing.
  out <- ifelse(x == "1", 1L, ifelse(x == "2", 0L, NA_integer_))
  out
}

items <- data.frame(
  smoked_before_13     = to_binary(get_col(lines, 145)),  # QN29
  smoked_daily_30d      = to_binary(get_col(lines, 150)),  # QN34
  drove_drinking        = to_binary(get_col(lines, 127)),  # QN11
  first_drink_before_13 = to_binary(get_col(lines, 156)),  # QN40
  binge_drink_30d       = to_binary(get_col(lines, 158)),  # QN42
  marijuana_before_13   = to_binary(get_col(lines, 161)),  # QN45
  cocaine_ever           = to_binary(get_col(lines, 164)),  # QN48
  glue_ever              = to_binary(get_col(lines, 166)),  # QN50
  meth_ever              = to_binary(get_col(lines, 168)),  # QN52
  ecstasy_ever           = to_binary(get_col(lines, 169)),  # QN53
  sex_before_13          = to_binary(get_col(lines, 174)),  # QN58
  sex_4plus_partners     = to_binary(get_col(lines, 175))   # QN59
)

grade_raw <- get_col(lines, 19)
grade <- factor(grade_raw, levels = c("1", "2", "3", "4"),
                labels = c("9", "10", "11", "12"))

sex_raw <- get_col(lines, 18)
sex <- factor(sex_raw, levels = c("1", "2"), labels = c("Female", "Male"))

weight  <- suppressWarnings(as.numeric(substr(lines, 358, 369)))
psu     <- trimws(substr(lines, 370, 374))
stratum <- trimws(substr(lines, 375, 378))

yrbs2005 <- data.frame(
  grade = grade,
  sex = sex,
  items,
  weight = weight,
  psu = psu,
  stratum = stratum
)

# Collins & Lanza (2010)'s analysis sample drops cases missing on `grade` as
# well as cases missing on every indicator: of the 13,917 raw records, 75
# have no grade recorded, and of the 3 cases missing on every indicator, 2
# of those do have a grade, giving 13917 - 75 - 2 = 13,840, matching the
# book's N exactly. Keep only that sample.
has_any_indicator <- Reduce(`|`, lapply(items, function(x) !is.na(x)))
yrbs2005 <- yrbs2005[!is.na(yrbs2005$grade) & has_any_indicator, ]

# Sanity check against Collins & Lanza (2010), Table 2.6 (marginal "Yes"
# proportions, unweighted, N = 13,840 responding).
expected <- c(smoked_before_13 = .15, smoked_daily_30d = .12,
              drove_drinking = .11, first_drink_before_13 = .26,
              binge_drink_30d = .25, marijuana_before_13 = .09,
              cocaine_ever = .08, glue_ever = .12, meth_ever = .06,
              ecstasy_ever = .06, sex_before_13 = .07,
              sex_4plus_partners = .17)
observed <- vapply(yrbs2005[names(items)], function(x) mean(x, na.rm = TRUE),
                    numeric(1))
stopifnot(nrow(yrbs2005) == 13840)
comparison <- round(rbind(observed, expected = expected[names(observed)]), 2)
print(comparison)
if (any(abs(comparison["observed", ] - comparison["expected", ]) > .01 + 1e-9))
  warning("A recoded item's marginal proportion doesn't match ",
          "Collins & Lanza (2010) Table 2.6 within rounding; check the ",
          "column position for that QN variable.")

usethis::use_data(yrbs2005, overwrite = TRUE)
