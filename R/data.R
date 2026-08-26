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

#' Simulated bystander intervention data (Liang & Park, Study 2 analogue)
#'
#' @description
#' **Synthetic data.** This dataset contains no real participant responses. It
#' is drawn from the estimated parameters of the published three-class latent
#' profile solution in Liang and Park's Study 2, so that the analysis can be
#' taught reproducibly. The authors' original data are not redistributed here.
#'
#' Five bystander-action indicators define three profiles: a *supportive-only*
#' profile (high emotional support, low everything else), a *disengaged*
#' profile (uniformly low), and a *broad responder* profile (uniformly high).
#'
#' @format A data frame with 300 rows and 19 variables:
#' \describe{
#'   \item{id}{Row identifier.}
#'   \item{age}{Age in years (18-70).}
#'   \item{male}{Sex, 0 = female, 1 = male. 3 values missing.}
#'   \item{org_intolerance_sh}{Organizational intolerance of sexual harassment
#'     (OITSH): perceived likelihood that the organization would respond to a
#'     harassment incident, 1-5, mean of 3 items.}
#'   \item{masc_job_context}{Masculine job-gender context: perceived percentage
#'     of males in the respondent's work unit and among their supervisors,
#'     0-1 slider.}
#'   \item{sh_experience}{Prior sexual-harassment experience, 1-4, mean of 28
#'     items. 1 missing.}
#'   \item{anger}{Anger at the harassment incident, 1-5, mean of 3 items.}
#'   \item{empathy}{Empathy for the target, 1-5, mean of 3 items.}
#'   \item{curb_expectancy}{Expectancy that sexual harassment can be curbed,
#'     1-5, mean of 4 items.}
#'   \item{confront}{LPA indicator: direct confrontation of the harasser, 1-5.}
#'   \item{distract}{LPA indicator: distraction or interruption, 1-5.}
#'   \item{support}{LPA indicator: emotional support to the target, 1-5.}
#'   \item{report}{LPA indicator: reporting to an authority, 1-5.}
#'   \item{discuss}{LPA indicator: speaking with the target afterwards, 1-5.}
#'   \item{harasser_aggression}{Distal outcome (BCH): the harasser's aggression
#'     toward the bystander after the intervention, single item, 1-5.}
#'   \item{target_gratitude}{Distal outcome (BCH): the target's gratitude
#'     toward the bystander after the intervention, single item, 1-5.
#'     1 missing.}
#'   \item{third_party_elevation}{Distal outcome (BCH): other third parties'
#'     moral elevation in response to the bystander's action, 1-5, mean of 4
#'     items. 1 missing.}
#'   \item{class_true}{Generating profile (1, 2, 3). Present only because the
#'     data are simulated; no real dataset would carry this.}
#'   \item{boundary}{\code{TRUE} for cases drawn from a blend of two profiles
#'     rather than one. These are the genuinely ambiguous respondents that give
#'     the data a realistic entropy; also a synthetic-only column.}
#' }
#'
#' @section Interpretation warning:
#' Fit statistics, standard errors, and p-values computed on this dataset
#' describe the simulation, not the original study. Do not cite them as
#' empirical findings about bystander intervention. Cite the paper for the
#' substantive results and this package for the simulated data.
#'
#' @section Known departures from the original:
#' A fraction of cases (flagged by \code{boundary}) are drawn from a blend of
#' two profiles, so the indicator block is deliberately not generated exactly
#' from the three-class model -- this is what gives the data a realistic
#' entropy instead of the artificially sharp separation you get by simulating
#' from a correctly specified model. Covariate-to-class logits are approximate
#' rather than exact; the indicators are independent of the covariates and
#' outcomes within class; item-level variables are not simulated. See
#' \code{vignette("liang_park_lpa")}.
#'
#' @source Simulated from the estimated three-class solution reported in
#'   Liang, Y., & Park, Y. (2025). A spectrum of bystander actions: Latent
#'   profile analysis of sexual harassment intervention behavior at work.
#'   \emph{Journal of Applied Psychology}. Advance online publication.
#'   \doi{10.1037/apl0001280}. The generator and the full generating-parameter
#'   tables are not distributed with the package.
#'
#' @examples
#' data(liang_park_sim)
#'
#' # The three profiles
#' aggregate(
#'   cbind(confront, distract, support, report, discuss) ~ class_true,
#'   data = liang_park_sim, FUN = mean
#' )
#'
#' \donttest{
#' # Blind recovery: free means, equal diagonal variances (as fitted originally)
#' ind <- c("confront", "distract", "support", "report", "discuss")
#' fit <- fit_mixture(liang_park_sim[, ind], n_classes = 3,
#'                    measurement = "continuous", n_init = 10)
#' table(class_assignments(fit), liang_park_sim$class_true)
#' }
"liang_park_sim"
