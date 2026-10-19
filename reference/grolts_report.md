# Answer the GRoLTS Reporting Checklist for a Trajectory Model

Walks the 21 items of the GRoLTS checklist (Guidelines for Reporting on
Latent Trajectory Studies; van de Schoot et al., 2017) for a fitted
growth model or a comparison of them, and sorts each item into one of
three groups: answered from the fit, which you must supply, and not
applicable to this model.

## Usage

``` r
grolts_report(object, ...)

# S3 method for class 'mixture_model'
grolts_report(object, ...)

# S3 method for class 'mixture_comparison'
grolts_report(object, k = NULL, ...)
```

## Arguments

- object:

  A
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  or
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  fit, or a comparison of them from
  [`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md).

- ...:

  Unused.

- k:

  For a comparison only: the number of classes of the final model.
  Default is the comparison's lowest-BIC model.

## Value

A data frame of class `grolts_report` with one row per item and columns
`item`, `question`, `status` (`"from the fit"`, `"you supply"` or
`"not applicable"`) and `answer`. Printing it groups the items by
status.

## Details

Van de Schoot et al. (2017) applied the checklist to 38 published
trajectory studies; only one reported the number of random starts (item
9), and none showed the mean trajectories of every model compared (item
14b). Most of what the checklist asks for is already on the fitted
object; this function gathers it, and names what the fit cannot answer:
the missing-data mechanism, the variables related to attrition, the
distribution of the observed variables, the plots, and which alternative
specifications were considered and why.

It deliberately prints no score. A count of items met invites treating
the items as interchangeable, which they are not.

A comparison from
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
answers the enumeration items (10, 11, 12) for every model; a single fit
answers them for itself only, and says so.

## References

van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., &
Vermunt, J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on
latent trajectory studies. *Structural Equation Modeling*, *24*(3),
451-467.
[doi:10.1080/10705511.2016.1247646](https://doi.org/10.1080/10705511.2016.1247646)

## See also

[`report_methods()`](https://pdvalencia.github.io/mixtureEM/reference/report_methods.md)
for the methods paragraph.

## Examples

``` r
set.seed(1)
y <- matrix(rnorm(400), 100) + rep(c(0, 1, 2, 3), each = 100)
fit <- fit_lcga(y, times = 4, n_classes = 2, family = "gaussian",
                n_init = 3, random_state = 1)
grolts_report(fit)
#> GRoLTS checklist (van de Schoot et al., 2017)
#> 
#> Answered from the fit
#> 1. Is the metric of time used in the statistical model reported?
#>     Time scores 0, 1, 2, 3 for occasions T1, T2, T3, T4.
#> 
#> 3c. Is a description provided of how missing data in the analyses were
#>     dealt with?
#>     No values were missing on the 4 occasions.
#> 
#> 5. Is the software mentioned?
#>     R package mixtureEM 0.6.0, R 4.6.1.
#> 
#> 9. Is information reported about the number of random start values and
#>     final iterations included?
#>     The search used 3 random starts; the best log-likelihood was
#>     reached by 3 of them. The algorithm converged in 105 iterations.
#> 
#> 12. Are the number of cases per class reported for each model?
#>     Final model only: 75.9%, 24.1%. Pass the comparison for every
#>     model.
#> 
#> 13. If classification of cases in a trajectory is the goal, is entropy
#>     reported?
#>     Relative entropy 0.353.
#> 
#> 15. Are characteristics of the final class solution numerically
#>     described?
#>     Growth-factor means by class are in summary(fit); class shares
#>     75.9%, 24.1%.
#> 
#> 16. Are the syntax files available?
#>     Models were estimated with the R package mixtureEM (version 0.6.0)
#>     in R 4.6.1. Random starts were drawn with random_state = 1. The
#>     model is refitted by: fit_lcga(indicator = y, n_classes = 2, times
#>     = 4, family = "gaussian", n_init = 3, random_state = 1)
#> 
#> 
#> You must supply
#> 2. Is information presented about the mean and variance of time within
#>     a wave?
#>     The model uses one fixed time score per occasion; report how much
#>     actual measurement times varied within each wave.
#> 
#> 3a. Is the missing data mechanism reported?
#>     The estimator assumes values are missing at random; argue whether
#>     that is plausible. mcar_test() tests the stronger MCAR assumption.
#> 
#> 3b. Is a description provided of what variables are related to
#>     attrition/missing data?
#>     Name the variables related to attrition or missingness.
#> 
#> 4. Is information about the distribution of the observed variables
#>     included?
#>     Describe the observed continuous outcome at each occasion (means,
#>     SDs, skew).
#> 
#> 6a. Are alternative specifications of within-class heterogeneity
#>     considered and clearly documented?
#>     This fit is a 2-class latent class growth analysis (linear growth,
#>     no within-class growth-factor variance). Report which other
#>     within-class specifications (e.g., LCGA against GMM) were fitted
#>     and why this one was kept.
#> 
#> 6b. Are alternative specifications of the between-class differences in
#>     variance-covariance matrix structure considered and clearly
#>     documented?
#>     This fit has no growth-factor covariance; its residual variances
#>     are free across classes. Report which other structures were
#>     considered.
#> 
#> 7. Are alternative shape/functional forms of the trajectories
#>     described?
#>     This fit uses linear growth. Report the other shapes fitted, if
#>     any.
#> 
#> 10. Are the model comparison (and selection) tools described from a
#>     statistical perspective?
#>     Pass the compare_longitudinal() result to answer this from the
#>     fits.
#> 
#> 11. Are the total number of fitted models reported, including a
#>     one-class solution?
#>     Pass the compare_longitudinal() result to answer this from the
#>     fits.
#> 
#> 14a. Is a plot included with the estimated mean trajectories of the
#>     final solution?
#>     plot(fit) draws the estimated mean trajectories; include it.
#> 
#> 14b. Are plots included with the estimated mean trajectories for each
#>     model?
#>     Plot each compared model's mean trajectories, not only the final
#>     one's.
#> 
#> 14c. Is a plot included of the estimated means of the final model with
#>     the observed individual trajectories, by class?
#>     Overlay the observed individual trajectories on each class's
#>     estimated mean.
#> 
#> 
#> Not applicable to this model
#> 8. If covariates have been used, can analyses still be replicated?
#>     No covariates were used.
#> 
```
