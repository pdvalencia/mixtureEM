
# mixtureEM

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdvalencia/mixtureEM/actions/workflows/R-CMD-check.yaml)

**mixtureEM** is latent class and latent profile analysis for applied
researchers: find the hidden subgroups in your data, describe them
(Collins & Lanza, 2010).

If you have ever run an LCA, assigned everyone to their most likely
class, and then run chi-squares on the labels, this package exists for
you: the step that is easiest to get wrong (relating classes to external
variables) is the one mixtureEM automates correctly, using the
bias-adjusted three-step estimators from the methodological literature
(Bakk et al., 2014; Vermunt, 2010).

And if you have never run an LCA, this package exists for you too: every
step ships with a sensible default and output that reads like the table
you would put in a paper, so you do not need the estimation theory
memorized to get correct answers out.

## Installation

``` r
# install.packages("pak")
pak::pkg_install("pdvalencia/mixtureEM")
```

## A complete analysis in a handful of calls

Using the bundled `ventura_leon` data — 16 yes/no infidelity items from
400 young adults (Ventura-León et al., 2025):

``` r
library(mixtureEM)

items <- ventura_leon[, 7:22]   # the 16 infidelity items
fit <- fit_mixture(items, n_classes = 4, n_init = 20, random_state = 1)
class_sizes(fit)
#>   class proportion n_expected n_modal
#> 1     1  0.4279122  171.16487     173
#> 2     2  0.2697511  107.90043     106
#> 3     3  0.1579645   63.18579      63
#> 4     4  0.1443723   57.74891      58
```

``` r
plot(fit, class_labels = c("Fidelity", "Affective interest",
                           "Infidelity", "Sexual desire"))
```

<img src="man/figures/README-example-plot-1.png" width="100%" />

Choosing how many classes to keep, reading the full item-response table,
relating classes to covariates with a bias correction, and describing
them by a distal outcome are each one function call away — the
[`ventura_leon` worked
example](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html)
walks through all of it on this exact fit.

## What it covers

| Analysis | Function | Worked example |
|----|----|----|
| Latent class analysis (binary / categorical / count items) | `fit_mixture()` | [`ventura_leon`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html) |
| Latent profile analysis (continuous items), missing data, multiple groups | `fit_mixture()` | [`janousch`](https://pdvalencia.github.io/mixtureEM/articles/janousch.html) |
| Choosing the number of classes (ICs, BLRT, fit diagnostics) | `compare_mixtures()`, `blrt()` | [`class_enumeration`](https://pdvalencia.github.io/mixtureEM/articles/class_enumeration.html) |
| Predictors of class membership (3-step ML) | `add_covariates()` | [`ventura_leon`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html) |
| Distal outcomes (BCH / ML) | `add_outcome()` | [`ventura_leon`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.html) |
| Repeated-measures LCA (trajectory classes) | `fit_rmlca()` | [`rmlca`](https://pdvalencia.github.io/mixtureEM/articles/rmlca.html) |
| Latent transition analysis, mover-stayer models | `fit_lta()` | [`lta`](https://pdvalencia.github.io/mixtureEM/articles/lta.html) |
| Latent class growth analysis, growth mixture models | `fit_lcga()`, `fit_gmm()` | [`growth_mixture`](https://pdvalencia.github.io/mixtureEM/articles/growth_mixture.html) |
| Complex survey designs (weights, strata, clusters) | every model | [`survey_lca`](https://pdvalencia.github.io/mixtureEM/articles/survey_lca.html) |

Also included: mixed measurement models (binary + continuous + count
indicators in one model), full-information missing-data handling
everywhere, absolute and local fit diagnostics (`absolute_fit()`,
`bivariate_residuals()`), classification diagnostics, and standard
errors that account for the classes being estimated rather than observed
(`?covariate_se`).

## Design principles

- **Conceptual workflow.** Choose a model, *then* examine covariates and
  outcomes on that fitted model — `add_covariates(fit, ...)` and
  `add_outcome(fit, ...)` never re-estimate the classes you already
  inspected.
- **Sound defaults.** Stepwise analyses apply the recommended bias
  correction and uncertainty-aware standard errors unless you say
  otherwise; where a default has a literature, the default follows it.
- **Readable output, reusable numbers.** Summaries print as formatted
  tables *and* return tidy data frames.
- **No dependency stack.** Base R only, so it installs anywhere R does.

## Citation

If you use mixtureEM in published research, please cite it as:

    Valencia, P. D. (2026). mixtureEM: Latent class and profile analysis via
    mixture modelling (Version 0.2.0) [R package].
    https://github.com/pdvalencia/mixtureEM

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

Ventura-León, J., Reyes, A., Valencia, P. D., Tocto-Muñoz, S.,
Gamboa-Melgar, G., Ruiz-Castro, J., & Lino-Cruz, C. (2025). Exploring
infidelity behavior patterns in a sample of Peruvian young adults: A
latent class analysis. *Journal of Marital and Family Therapy*, *51*(4),
e70066. <https://doi.org/10.1111/jmft.70066>

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
<https://doi.org/10.1093/pan/mpq025>
