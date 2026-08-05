# Prepare the bundled `janousch` dataset from Janousch et al. (2022),
# "Resilience profiles across context: A latent profile analysis in a
# German, Greek, and Swiss sample of adolescents", PLOS ONE 17(1): e0263089.
# https://doi.org/10.1371/journal.pone.0263089
#
# Source: figshare supplementary file, deposited by the paper's authors
# under a CC BY 4.0 license (confirmed via the figshare API:
# https://api.figshare.com/v2/articles/19081099).
# https://figshare.com/articles/dataset/Dataset_/19081099
#
# Raw file has 13 columns: country, gender, age, 5 Resilience Scale for
# Adolescents (READ) subscale means (pc, sc, ss, sr, fc), 3 Hopkins Symptom
# Checklist scores (hscl total, anx, dep), and migration background (mig).
# There is substantial item-level missingness (see paper sec. 2.4: analyses
# used Mplus's MLR estimator, i.e. full-information ML over the available
# indicators, not listwise deletion) -- kept as NA here so it can be passed
# to `fit_mixture()`'s `"continuous_nan"` measurement type, which handles
# missing continuous indicators the same way.

raw <- haven::read_sav("data-raw/janousch/janousch_raw.sav")

country <- factor(haven::as_factor(raw$country),
                   levels = c("Switzerland", "Germany", "Greece"))
gender <- factor(haven::as_factor(raw$gender),
                  levels = c("male", "female", "other"),
                  labels = c("Male", "Female", "Other"))

mig_num <- as.numeric(raw$mig)
migration_background <- factor(ifelse(mig_num < 0, NA, mig_num),
                                levels = c(1, 0),
                                labels = c("Native", "Migration background"))

janousch <- data.frame(
  country = country,
  gender = gender,
  age = as.numeric(raw$age),
  migration_background = migration_background,
  anxiety = as.numeric(raw$anx),
  depression = as.numeric(raw$dep),
  personal_competence = as.numeric(raw$pc),
  social_competence = as.numeric(raw$sc),
  structured_style = as.numeric(raw$ss),
  social_resources = as.numeric(raw$sr),
  family_cohesion = as.numeric(raw$fc)
)

# Sanity checks against Table 1 (per-country subscale means) and the
# reported sample sizes/composition (sec. 2.1) of Janousch et al. (2022).
stopifnot(nrow(janousch) == 1160)
stopifnot(all(table(janousch$country) ==
                c(Switzerland = 375, Germany = 346, Greece = 439)))
stopifnot(all(table(janousch$gender) ==
                c(Male = 565, Female = 565, Other = 2)))

country_means <- function(var) {
  round(tapply(janousch[[var]], janousch$country, mean, na.rm = TRUE), 2)
}
published <- list(
  personal_competence = c(Switzerland = 3.89, Germany = 3.88, Greece = 3.93),
  social_competence   = c(Switzerland = 4.02, Germany = 3.95, Greece = 4.12),
  structured_style    = c(Switzerland = 3.67, Germany = 3.61, Greece = 3.75),
  social_resources    = c(Switzerland = 4.50, Germany = 4.50, Greece = 4.57),
  family_cohesion     = c(Switzerland = 4.33, Germany = 4.32, Greece = 4.26),
  anxiety             = c(Switzerland = 1.95, Germany = 1.89, Greece = 1.77),
  depression          = c(Switzerland = 1.78, Germany = 1.80, Greece = 1.71)
)
for (var in names(published)) {
  stopifnot(all(abs(country_means(var)[names(published[[var]])] -
                       published[[var]]) < 0.01))
}

usethis::use_data(janousch, overwrite = TRUE)
