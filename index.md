# mixtureEM

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml)

Latent class and latent profile analysis for applied researchers.

You already suspect your data holds distinct subgroups — different ways
people respond to a workplace survey, different patterns of health-risk
behavior, different trajectories through an intervention. You know these
hidden subgroups exist; you just need a principled way to coax them out
of the data rather than guessing at cutpoints. **mixtureEM** does that
(Collins & Lanza, 2010), and describes the groups it finds in plain,
reusable output.

- If you have ever run an LCA, assigned everyone to their most likely
  class, and then run ANOVAs or chi-squares on the labels, this package
  exists for you. The step that is easiest to get wrong — relating
  classes to external variables — is the one mixtureEM automates
  correctly, using the bias-adjusted three-step estimators from the
  methodological literature (Bakk et al., 2014; Vermunt, 2010).
- If you have never run an LCA, this package exists for you too. Every
  step ships with a sensible default, and the output reads like the
  table you would put in a paper, so you do not need the estimation
  theory memorized to get correct answers out.

------------------------------------------------------------------------

## Installation

You can install the development version from GitHub using the
[pak](https://pak.r-lib.org/) package:

``` r

# install.packages("pak")
pak::pkg_install("pdvalencia/mixtureEM")
```

------------------------------------------------------------------------

## Quick start

A complete analysis takes three steps: decide how many classes to keep,
fit the model and describe the classes it finds, and relate those
classes to the variables you actually care about. All three use the
bundled `ventura_leon` data — 16 yes/no infidelity items from 400 young
adults (run
[`?ventura_leon`](https://pdvalencia.github.io/mixtureEM/reference/ventura_leon.md)
for dataset details).

**Step 1: How many classes fit best?**

[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
fits a range of class counts and reports the information criteria and
diagnostics side by side. At the package’s default twenty random
restarts per class count, six models is real work, so it is shown here
rather than run on every render of this page — the [`class_enumeration`
vignette](https://pdvalencia.github.io/mixtureEM/articles/class_enumeration.html)
walks through reading its output on this exact data:

``` r

library(mixtureEM)
data(ventura_leon)
compare_mixtures(ventura_leon[, 7:22], measurement = "binary", k_range = 1:6)
```

**Step 2: Fit the model and describe the classes.**

Four classes turn out to be the interpretable choice here:

``` r

library(mixtureEM)

items <- ventura_leon[, 7:22]   # the 16 infidelity items
fit <- fit_mixture(items, n_classes = 4, measurement = "binary",
                   n_init = 20, random_state = 1)
class_sizes(fit)
#>   class proportion n_expected n_modal
#> 1     1  0.4276388  171.05551     173
#> 2     2  0.2699044  107.96177     106
#> 3     3  0.1581894   63.27577      63
#> 4     4  0.1442674   57.70695      58
```

``` r

plot(fit, class_labels = c("Fidelity", "Affective interest",
                           "Infidelity", "Sexual desire"))
```

![Item-response probabilities for the four fitted classes across the 16
infidelity items, one line per
class.](reference/figures/README-example-plot-1.png)

**Step 3: Relate the classes to covariates.**

Who ends up in which class?
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
applies the bias-adjusted three-step correction (Bakk et al., 2014;
Vermunt, 2010) automatically, so the classes you just interpreted do not
shift under you:

``` r

fit_cov <- add_covariates(fit, ~ sex + age, data = ventura_leon)
summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                               OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                1.178  [    0.234,     5.924]     0.842
#>   sex.Female               2.641  [    1.104,     6.318]     0.029
#>   age                      0.942  [    0.884,     1.003]     0.063
#> 
#> Class 3 ON
#>   Intercept                0.315  [    0.104,     0.953]     0.041
#>   sex.Female               0.342  [    0.183,     0.641]    < .001
#>   age                      1.033  [    0.997,     1.071]     0.074
#> 
#> Class 4 ON
#>   Intercept                0.247  [    0.078,     0.779]     0.017
#>   sex.Female               0.656  [    0.331,     1.299]     0.226
#>   age                      1.024  [    0.985,     1.065]     0.227
#> 
#> OMNIBUS TEST PER COVARIATE (effect across all classes)
#> ---------------------------------------------------------
#>                          Wald Chi2   df  P-Value
#>   sex                       23.066    3    < .001
#>   age                        9.645    3     0.022
#>   Note: a non-significant test beside large coefficients can be the
#>         Hauck-Donner effect; confirm with wald_omnibus_test().
#> =========================================================
```

> Reading the full item-response table and describing classes by a
> distal outcome are each one function call away too. The
> [`ventura_leon` worked
> example](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html)
> walks through all of it on this exact fit.

------------------------------------------------------------------------

## What it covers

| Analysis | Function | Worked example |
|----|----|----|
| Latent class analysis (binary / categorical / count items) | [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md) | [`ventura_leon`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html) |
| Latent profile analysis (continuous items), class enumeration, covariates and distal outcomes | [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md) | [`liang_park_lpa`](https://pdvalencia.github.io/mixtureEM/articles/liang_park_lpa.html) |
| Multiple-group LCA: do the classes mean the same thing in every group before comparing their prevalences? | `fit_mixture(group = )` | [`mglca_yrbs`](https://pdvalencia.github.io/mixtureEM/articles/mglca_yrbs.html) |
| Choosing the number of classes (ICs, BLRT, fit diagnostics) | [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md), [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md) | [`class_enumeration`](https://pdvalencia.github.io/mixtureEM/articles/class_enumeration.html) |
| Predictors of class membership (3-step ML) | [`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md) | [`ventura_leon`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html) |
| Distal outcomes (BCH / ML) | [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md) | [`liang_park_lpa`](https://pdvalencia.github.io/mixtureEM/articles/liang_park_lpa.html) |
| Repeated-measures LCA: trajectory classes over time, e.g. “Stable Low”, “Escalating”, “Persistent” | [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md) | [`rmlca`](https://pdvalencia.github.io/mixtureEM/articles/rmlca.html) |
| Latent transition analysis, mover-stayer models: if someone is Depressed at Wave 1, what is the probability they are Not Depressed at Wave 2? | [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md) | [`lta`](https://pdvalencia.github.io/mixtureEM/articles/lta.html) |
| Latent class growth analysis, growth mixture models | [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md), [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md) | [`growth_mixture`](https://pdvalencia.github.io/mixtureEM/articles/growth_mixture.html) |
| Complex survey designs (weights, strata, clusters) | every model | [`survey_lca`](https://pdvalencia.github.io/mixtureEM/articles/survey_lca.html) |

Also included: mixed measurement models (binary + continuous + count
indicators in one model), full-information missing-data handling
everywhere, absolute and local fit diagnostics
([`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)),
classification diagnostics, and standard errors that account for the
classes being estimated rather than observed
([`?covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md)).

------------------------------------------------------------------------

## Design principles

- **Conceptual workflow.** Choose a model, *then* examine covariates and
  outcomes on that fitted model — `add_covariates(fit, ...)` and
  `add_outcome(fit, ...)` never re-estimate the classes you already
  inspected.
- **Guards against the classify-and-analyze trap.** The most common
  mistake in mixture modeling — one you may already have made without
  realizing it — is assigning each case to its most likely class and
  then running an ordinary ANOVA or chi-square on the labels. That
  ignores classification uncertainty and biases results toward the null.
  Stepwise analyses apply the recommended bias correction and
  uncertainty-aware standard errors unless you say otherwise; where a
  default has a literature, the default follows it.
- **Takes your survey design seriously.** Large surveys are rarely a
  simple random sample. When your data carries sampling weights, strata,
  and clusters, mixtureEM accounts for them — it will not just pretend
  your complex survey is a simple random sample.
- **Readable output, reusable numbers.** Summaries print as formatted
  tables *and* return tidy data frames, not deeply nested lists of
  matrices you have to pick apart by hand.
- **No dependency stack.** Built using base R only, meaning it installs
  cleanly anywhere R does.

------------------------------------------------------------------------

## Citation

If you use mixtureEM in published research, please cite it as:

> Valencia, P. D. (2026). mixtureEM: Latent class and profile analysis
> via mixture modelling (Version 0.3.0) \[R package\].
> <https://github.com/pdvalencia/mixtureEM>

## Acknowledgements

**mixtureEM** draws strong inspiration from the Python package
[**StepMix**](https://github.com/Labo-Lacourse/stepmix), which pioneered
open-source, bias-adjusted multi-step estimation of generalised mixture
models with external variables. If you make use of these methods, please
also consider citing the StepMix paper:

> Morin, S., Legault, R., Laliberte, F., Bakk, Z., Giguère, C.-E., de la
> Sablonnière, R., & Lacourse, E. (2025). StepMix: A Python package for
> pseudo-likelihood estimation of generalised mixture models with
> external variables. *Journal of Statistical Software*, *113*(8), 1–39.
> <https://doi.org/10.18637/jss.v113.i08>

## References

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, *22*(4), 520–540.
<https://doi.org/10.1093/pan/mpu003>

Collins, L. M., & Lanza, S. T. (2010). *Latent class and latent
transition analysis: With applications in the social, behavioral, and
health sciences*. Wiley.

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
<https://doi.org/10.1093/pan/mpq025>
