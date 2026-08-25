# Simulated bystander intervention data (Liang & Ark, Study 2 analogue)

**Synthetic data.** This dataset contains no real participant responses.
It is drawn from the estimated parameters of the published three-class
latent profile solution in Liang and Ark's Study 2, so that the analysis
can be taught reproducibly. The authors' original data are not
redistributed here.

Five bystander-action indicators define three profiles: a
*supportive-only* profile (high emotional support, low everything else),
a *disengaged* profile (uniformly low), and a *broad responder* profile
(uniformly high).

## Usage

``` r
liang_ark_sim
```

## Format

A data frame with 300 rows and 19 variables:

- id:

  Row identifier.

- age:

  Age in years (18-70).

- male:

  Sex, 0 = female, 1 = male. 3 values missing.

- org_intolerance_sh:

  Organizational intolerance of sexual harassment (OITSH): perceived
  likelihood that the organization would respond to a harassment
  incident, 1-5, mean of 3 items.

- masc_job_context:

  Masculine job-gender context: perceived percentage of males in the
  respondent's work unit and among their supervisors, 0-1 slider.

- sh_experience:

  Prior sexual-harassment experience, 1-4, mean of 28 items. 1 missing.

- anger:

  Anger at the harassment incident, 1-5, mean of 3 items.

- empathy:

  Empathy for the target, 1-5, mean of 3 items.

- curb_expectancy:

  Expectancy that sexual harassment can be curbed, 1-5, mean of 4 items.

- confront:

  LPA indicator: direct confrontation of the harasser, 1-5.

- distract:

  LPA indicator: distraction or interruption, 1-5.

- support:

  LPA indicator: emotional support to the target, 1-5.

- report:

  LPA indicator: reporting to an authority, 1-5.

- discuss:

  LPA indicator: speaking with the target afterwards, 1-5.

- harasser_aggression:

  Distal outcome (BCH): the harasser's aggression toward the bystander
  after the intervention, single item, 1-5.

- target_gratitude:

  Distal outcome (BCH): the target's gratitude toward the bystander
  after the intervention, single item, 1-5. 1 missing.

- third_party_elevation:

  Distal outcome (BCH): other third parties' moral elevation in response
  to the bystander's action, 1-5, mean of 4 items. 1 missing.

- class_true:

  Generating profile (1, 2, 3). Present only because the data are
  simulated; no real dataset would carry this.

- boundary:

  `TRUE` for cases drawn from a blend of two profiles rather than one.
  These are the genuinely ambiguous respondents that give the data a
  realistic entropy; also a synthetic-only column.

## Source

Simulated from the estimated three-class solution reported in Liang &
Ark, *A spectrum of bystander actions: Latent profile analysis of sexual
harassment intervention*. The generator and the full
generating-parameter tables are not distributed with the package.

## Interpretation warning

Fit statistics, standard errors, and p-values computed on this dataset
describe the simulation, not the original study. Do not cite them as
empirical findings about bystander intervention. Cite the paper for the
substantive results and this package for the simulated data.

## Known departures from the original

A fraction of cases (flagged by `boundary`) are drawn from a blend of
two profiles, so the indicator block is deliberately not generated
exactly from the three-class model – this is what gives the data a
realistic entropy instead of the artificially sharp separation you get
by simulating from a correctly specified model. Covariate-to-class
logits are approximate rather than exact; the indicators are independent
of the covariates and outcomes within class; item-level variables are
not simulated. See
[`vignette("liang_ark_lpa")`](https://pdvalencia.github.io/mixtureEM/articles/liang_ark_lpa.md).

## Examples

``` r
data(liang_ark_sim)

# The three profiles
aggregate(
  cbind(confront, distract, support, report, discuss) ~ class_true,
  data = liang_ark_sim, FUN = mean
)
#>   class_true confront distract  support   report  discuss
#> 1          1 1.756522 2.284783 4.208696 1.933333 2.430435
#> 2          2 1.474638 1.758152 1.877717 1.380435 1.521739
#> 3          3 3.555556 3.577957 4.599462 3.759857 3.500000

# \donttest{
# Blind recovery: free means, equal diagonal variances (as fitted originally)
ind <- c("confront", "distract", "support", "report", "discuss")
fit <- fit_mixture(liang_ark_sim[, ind], n_classes = 3,
                   measurement = "continuous", n_init = 10)
table(class_assignments(fit), liang_ark_sim$class_true)
#>    
#>       1   2   3
#>   1 102   8  17
#>   2   7  80   2
#>   3   6   4  74
# }
```
