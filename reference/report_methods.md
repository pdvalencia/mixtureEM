# Write the Methods Paragraph for a Fitted Mixture Model

Returns the paragraphs a methods section needs to describe a fitted
model: the sample, the estimator and its random starts, missing data,
the fit indices, classification quality, any structural model, and the
line of R that refits it. Every number is read off the fitted object, so
the text is true of the fit it was written from and of no other.

## Usage

``` r
report_methods(object, ...)

# S3 method for class 'mixture_model'
report_methods(object, ...)

# S3 method for class 'lta_model'
report_methods(object, ...)

# S3 method for class 'mixture_comparison'
report_methods(object, k = NULL, ...)
```

## Arguments

- object:

  A fitted model
  ([`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md),
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md),
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md))
  or a comparison from
  [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
  or
  [`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md).

- ...:

  Unused.

- k:

  For a comparison only: the number of classes of the model the
  paragraphs describe. Default is the comparison's lowest-BIC model.

## Value

A character vector of class `methods_report`, one paragraph per element,
with the cited works in `attr(, "references")`. Printing it wraps the
paragraphs and lists the references.

## Details

**Why this exists.** Sorgente et al. (2025) coded 313 latent class
analyses published in psychology between 2014 and 2024 against the
published tutorials: 23% reported the estimator, 16% the random starts,
17% how missing data were handled, and 1.6% their analysis code. Every
one of those items is a value the package already stores; this function
puts them in the order and the wording a methods section uses.

**What it does not do.** It reports and never grades. No sentence calls
an entropy, a class size or a replication count adequate: where the
literature attaches a cutoff, the cutoff is discussed in the help page
that computes the number (for entropy,
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)).
It does not say why the number of classes is substantively right, what
the classes mean, or whether the sample suits the question; those are
the author's.

**Passing a comparison** (from
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
or
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md))
adds the enumeration paragraph: the range of K, whether the one-class
model was in it, the indices reported, and which K had the lowest BIC.
The other paragraphs then describe the model with `k` classes, by
default the lowest-BIC one.

Citations are returned as the `references` attribute and printed once at
the end, so the paragraphs serve both in-text and numbered citation
styles.

## References

Sorgente, A., Caliciuri, R., Robba, M., Lanz, M., & Zumbo, B. D. (2025).
A systematic review of latent class analysis in psychology: Examining
the gap between guidelines and research practice. *Behavior Research
Methods*, *57*, 301.

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

Weller, B. E., Bowen, N. K., & Faubert, S. J. (2020). Latent class
analysis: A guide to best practice. *Journal of Black Psychology*,
*46*(4), 287-311.
[doi:10.1177/0095798420930932](https://doi.org/10.1177/0095798420930932)

Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions
about latent class analysis. *Translational Issues in Psychological
Science*, *4*(4), 440-461.
[doi:10.1037/tps0000176](https://doi.org/10.1037/tps0000176)

## See also

[`grolts_report()`](https://pdvalencia.github.io/mixtureEM/reference/grolts_report.md)
for the trajectory-model checklist.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 5,
                   random_state = 1)
report_methods(fit)
#> The analysis sample comprised 100 cases. Classes were defined by 6
#> binary indicators.
#> 
#> A 2-class latent class model was estimated by maximum likelihood with
#> the EM algorithm (Dempster, Laird, & Rubin, 1977), with weak priors
#> that keep estimates off the boundary of the parameter space (Galindo
#> Garre & Vermunt, 2006): 1 pseudo-observation for the class proportions
#> and 1 pseudo-observation for the item-response probabilities. The
#> search used 5 random starts; the best log-likelihood was reached by 5
#> of them. The algorithm converged in 180 iterations.
#> 
#> The model's log-likelihood was -405.31 with 13 free parameters (AIC =
#> 836.62, BIC = 870.49, sample-size-adjusted BIC = 829.43).
#> 
#> Relative entropy was 0.427. The smallest class held 34.7% of the sample
#> (about 35 cases). The average posterior probability of the assigned
#> class ranged from 0.758 to 0.839 across classes. These describe how
#> cleanly cases are classified; they are not criteria for the number of
#> classes (Masyn, 2013).
#> 
#> Models were estimated with the R package mixtureEM (version 0.6.0) in R
#> 4.6.1. Random starts were drawn with random_state = 1. The model is
#> refitted by: fit_mixture(indicators = X, n_classes = 2, measurement =
#> "binary", n_init = 5, random_state = 1)
#> 
#> References
#> Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum
#>   likelihood from incomplete data via the EM algorithm. Journal of the
#>   Royal Statistical Society: Series B, 39(1), 1-38.
#> Galindo Garre, F., & Vermunt, J. K. (2006). Avoiding boundary estimates
#>   in latent class analysis by Bayesian posterior mode estimation.
#>   Behaviormetrika, 33(1), 43-59. https://doi.org/10.2333/bhmk.33.43
#> Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
#>   In T. D. Little (Ed.), The Oxford Handbook of Quantitative Methods
#>   (Vol. 2, pp. 551-611). Oxford University Press.
```
