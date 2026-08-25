# ==============================================================================
# Bundled example datasets
# ==============================================================================

#' Youth Risk Behavior Survey, 2005 (national public-use file)
#'
#' Twelve dichotomous health-risk-behavior items for U.S. high school
#' students, reproducing the empirical illustration used throughout Collins &
#' Lanza (2010), \emph{Latent Class and Latent Transition Analysis}
#' (sections 1.8.1.2, 2.4, 4.6.4 and 5.12.2). The items are CDC's own
#' precomputed dichotomous recodes (the `QN*` variables in the CDC file), so
#' they match CDC's and Collins & Lanza's coding exactly. `grade` is included
#' to reproduce the multiple-group analysis in sec. 5.12.2; `weight`, `psu`
#' and `stratum` reproduce the complex sample design and can be passed to
#' [`fit_mixture()`]'s `weights=`, `cluster=` and `strata=` arguments.
#'
#' @format A data frame with 13,840 rows (respondents) and 17 columns. This is
#'   Collins & Lanza's analysis sample: cases missing on `grade`, and the
#'   cases missing on every indicator, are dropped from the raw CDC file's
#'   13,917 records.
#' \describe{
#'   \item{grade}{Factor, grade in school: `"9"`, `"10"`, `"11"`, `"12"`.}
#'   \item{sex}{Factor: `"Female"`, `"Male"`.}
#'   \item{smoked_before_13}{Smoked a whole cigarette before age 13 (0/1).}
#'   \item{smoked_daily_30d}{Smoked cigarettes daily for 30 days (0/1).}
#'   \item{drove_drinking}{Drove after drinking alcohol, past 30 days (0/1).}
#'   \item{first_drink_before_13}{Had a first drink of alcohol before age 13
#'     (0/1).}
#'   \item{binge_drink_30d}{Had 5 or more drinks in a row, past 30 days
#'     (0/1).}
#'   \item{marijuana_before_13}{Tried marijuana before age 13 (0/1).}
#'   \item{cocaine_ever}{Used cocaine, ever (0/1).}
#'   \item{glue_ever}{Sniffed glue or inhalants, ever (0/1).}
#'   \item{meth_ever}{Used methamphetamines, ever (0/1).}
#'   \item{ecstasy_ever}{Used ecstasy, ever (0/1).}
#'   \item{sex_before_13}{Had sexual intercourse before age 13 (0/1).}
#'   \item{sex_4plus_partners}{Had sexual intercourse with four or more
#'     people, ever (0/1).}
#'   \item{weight}{Analysis (sampling) weight.}
#'   \item{psu}{Primary sampling unit identifier.}
#'   \item{stratum}{Sampling stratum identifier.}
#' }
#'
#' @source Centers for Disease Control and Prevention (2006). Youth Risk
#'   Behavior Survey, 2005 national public-use data set. Downloaded from
#'   \url{https://ftp.cdc.gov/pub/data/yrbs/2005/}. YRBS national data are
#'   produced by a U.S. federal agency and are in the public domain (17
#'   U.S.C. sec. 105); no permission is required to use or redistribute them.
#'   See `data-raw/yrbs2005.R` for the exact recoding script.
#'
#' @references Collins, L. M., & Lanza, S. T. (2010). \emph{Latent Class and
#'   Latent Transition Analysis: With Applications in the Social,
#'   Behavioral, and Health Sciences}. Wiley. Table 2.6 reports the marginal
#'   "Yes" proportions this recoding was checked against.
"yrbs2005"

#' Infidelity behavior patterns in Peruvian young adults
#'
#' Sixteen binary indicators from the Multidimensional Infidelity
#' Inventory-S (MII-S; Romero-Palencia et al. 2007) plus six sociodemographic
#' variables, from a latent class analysis of infidelity-related thoughts and
#' behaviors in a sample of 400 Peruvian young adults. Reproduces the
#' four-class solution (infidelity, sexual desire, affective interest,
#' fidelity) of Ventura-Leon et al. (2025). The original paper related
#' classes to covariates by assigning each respondent to their most likely
#' class and then running chi-square/ANOVA tests on the observed class
#' labels -- a classify-and-analyze approach that ignores classification
#' uncertainty. `fit_mixture()`'s `predictors=` argument with `n_steps = 3`
#' and ML bias correction is used instead in the package vignette, which
#' accounts for that uncertainty when testing the same covariates (sex, age,
#' sexual orientation, relationship duration).
#'
#' @format A data frame with 400 rows (respondents) and 22 columns:
#' \describe{
#'   \item{sex}{Factor: `"Male"`, `"Female"`.}
#'   \item{age}{Age in years.}
#'   \item{birthplace}{Factor: `"Lima"`, `"Outside Lima"`.}
#'   \item{sexual_orientation}{Factor: `"Heterosexual"`, `"Not
#'     heterosexual"`.}
#'   \item{relationship_duration}{Factor: `"Short"`, `"Long"`, dichotomized
#'     at the sample median (50 months), per the original study.}
#'   \item{relationship_type}{Factor: `"Married"`, `"Cohabiting"`, `"In
#'     love"`, `"Engaged"`, `"Dating"`.}
#'   \item{flirting}{"I have flirted with other people besides my partner."
#'     (0/1)}
#'   \item{romantic_partners}{"I have had another romantic partner(s)."
#'     (0/1)}
#'   \item{emotional_bond}{"I have formed an emotional bond with someone
#'     else besides my partner." (0/1)}
#'   \item{romantic_involvement}{"I have been romantically involved with
#'     another person/other people." (0/1)}
#'   \item{loved_another}{"I have loved another person(s) besides my
#'     partner." (0/1)}
#'   \item{in_love}{"I have fallen in love with someone else besides my
#'     partner." (0/1)}
#'   \item{thoughts}{"I have thought about someone else besides my partner."
#'     (0/1)}
#'   \item{interest}{"I have been interested in another person(s) besides my
#'     partner." (0/1)}
#'   \item{sexual_relations}{"I have had sexual relations with another
#'     person/other people besides my partner." (0/1)}
#'   \item{sexual_contact}{"I have had sexual contact with someone else
#'     besides my partner." (0/1)}
#'   \item{desired_relations}{"I have desired to have sexual relations with
#'     another person(s) besides my partner." (0/1)}
#'   \item{desired_contact}{"I have desired to have sexual contact with
#'     another person(s) besides my partner." (0/1)}
#'   \item{sexual_fantasies}{"I have wanted to fulfill my sexual fantasies
#'     with someone else besides my partner." (0/1)}
#'   \item{attraction}{"I have felt attracted to another person(s) besides
#'     my partner." (0/1)}
#'   \item{had_sex}{"I have had sex with another person/other people besides
#'     my partner." (0/1)}
#'   \item{desired_sex}{"I have desired to have sex with another person(s)
#'     besides my partner." (0/1)}
#' }
#'
#' @source OSF repository \url{https://osf.io/8csr9/}. The OSF project has no
#'   license flag set; Pablo D. Valencia, a co-author of the source paper and
#'   of this package, has authorized bundling this dataset in mixtureEM on
#'   behalf of the author team. See `data-raw/ventura_leon.R` for the exact
#'   recoding script.
#'
#' @references Ventura-Leon, J., Reyes, A., Valencia, P. D., Tocto-Munoz, S.,
#'   Gamboa-Melgar, G., Ruiz-Castro, J., & Lino-Cruz, C. (2025). Exploring
#'   infidelity behavior patterns in a sample of Peruvian young adults: A
#'   latent class analysis. \emph{Journal of Marital and Family Therapy},
#'   \emph{51}, e70066. \doi{10.1111/jmft.70066}. Tables 1-4 report the
#'   sociodemographic, item-response and covariate numbers this recoding and
#'   the package vignette are checked against.
"ventura_leon"

#' Resilience profiles in German, Greek, and Swiss adolescents
#'
#' Seven continuous indicators -- anxiety and depression (Hopkins Symptom
#' Checklist) and five protective-factor subscale means from the Resilience
#' Scale for Adolescents (personal competence, social competence, structured
#' style, social resources, family cohesion) -- for 1,160 seventh-graders
#' from Germany, Greece, and Switzerland, together with gender and migration
#' background. Reproduces the latent profile analysis of Janousch et al.
#' (2022), who fit a *separate* LPA per country (finding 3 profiles for
#' Switzerland and 4 for Germany and Greece) and then tested whether the
#' profiles were measurement-invariant across countries. Item-level
#' indicators have substantial missingness (see `@format`); the original
#' study used full-information ML (Mplus MLR) rather than listwise deletion,
#' which corresponds to `fit_mixture()`'s `"continuous_nan"` measurement
#' type.
#'
#' @format A data frame with 1,160 rows (students) and 11 columns:
#' \describe{
#'   \item{country}{Factor: `"Switzerland"`, `"Germany"`, `"Greece"`.}
#'   \item{gender}{Factor: `"Male"`, `"Female"`, `"Other"` (2 cases).}
#'   \item{age}{Age in years (76 missing).}
#'   \item{migration_background}{Factor: `"Native"`, `"Migration
#'     background"` (21 missing).}
#'   \item{anxiety}{Hopkins Symptom Checklist anxiety subscale mean, 1-4
#'     (143 missing).}
#'   \item{depression}{Hopkins Symptom Checklist depression subscale mean,
#'     1-4 (177 missing).}
#'   \item{personal_competence,social_competence,structured_style,
#'     social_resources,family_cohesion}{Resilience Scale for Adolescents
#'     (READ) subscale means, 1-5 (84-132 missing per subscale).}
#' }
#'
#' @source figshare supplementary file
#'   \url{https://figshare.com/articles/dataset/Dataset_/19081099}, deposited
#'   by the paper's authors under a CC BY 4.0 license. See
#'   `data-raw/janousch.R` for the exact recoding script.
#'
#' @references Janousch, C., Anyan, F., Kassis, W., Morote, R., Hjemdal, O.,
#'   Sidler, P., Graf, U., Rietz, C., Chouvati, R., & Govaris, C. (2022).
#'   Resilience profiles across context: A latent profile analysis in a
#'   German, Greek, and Swiss sample of adolescents. \emph{PLOS ONE},
#'   \emph{17}(1), e0263089. \doi{10.1371/journal.pone.0263089}. Table 1
#'   reports the per-country subscale means this recoding is checked
#'   against; Table 2 and Figs 1-3 report the profile solutions and Table 3
#'   the profile means the package vignette is checked against.
"janousch"
