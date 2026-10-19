# Reporting a Mixture Model

``` r

library(mixtureEM)
```

## Why a reporting function

Sorgente et al. (2025) read 313 latent class analyses published in
psychology between 2014 and 2024 and checked them against what the
tutorials say a paper should report. Most papers left out the basics:
23% named the estimator, 16% said how many random starts they used, 17%
said how missing data were handled, and 1.6% shared their code. Van de
Schoot et al. (2017) found the same pattern in trajectory studies: of 38
papers, one reported its number of random starts.

None of this is hard to report. The fitted model already holds every one
of those numbers; they are just spread over
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md)
and the comparison table.
[`report_methods()`](https://pdvalencia.github.io/mixtureEM/reference/report_methods.md)
collects them into the paragraphs a methods section needs, and
[`grolts_report()`](https://pdvalencia.github.io/mixtureEM/reference/grolts_report.md)
walks the GRoLTS checklist for trajectory models. Neither function
computes anything new, and neither grades the result: no sentence calls
an entropy “good” or a class “too small”. Those judgements need the
substantive context the function does not have.

## A latent class analysis, from enumeration to methods section

The `ventura_leon` data hold 16 binary infidelity items for 400 young
adults (see
[`vignette("ventura_leon")`](https://pdvalencia.github.io/mixtureEM/articles/ventura_leon.md)).
Compare one to four classes:

``` r

items <- as.matrix(ventura_leon[, 7:22])
cmp <- compare_mixtures(items, k_range = 1:4, measurement = "binary",
                        n_init = 20, random_state = 1)
#> Running Model Selection across K = 1 to 4...
#> 
#> Fitting 1-class model...
#> Fitting 2-class model...
#> Fitting 3-class model...
#> Fitting 4-class model...
#> 
#> === Model Selection Summary ===
#>  Classes         LL Params      BIC    SABIC     AIC3 Entropy
#>        1 -3905.807      16 7907.478 7856.709 7859.615   1.000
#>        2 -2900.078      33 5997.875 5893.163 5899.156   0.948
#>        3 -2654.569      50 5608.711 5450.057 5459.138   0.934
#>        4 -2531.141      67 5463.710 5251.115 5263.282   0.907
#> 
#> -> Best model according to BIC: 4 classes
#>    print(x, columns = "all") adds AIC, CAIC, AWE, ICL, BF and cmP.
```

Passing the comparison to
[`report_methods()`](https://pdvalencia.github.io/mixtureEM/reference/report_methods.md)
writes up the model with the lowest BIC, plus a paragraph on the
enumeration itself: the range of K, whether the one-class model was in
it (a common omission), the indices reported and which K the BIC picked.

``` r

report_methods(cmp)
#> The analysis sample comprised 400 cases. Classes were defined by 16
#> binary indicators.
#> 
#> A 4-class latent class model was estimated by maximum likelihood with
#> the EM algorithm (Dempster, Laird, & Rubin, 1977), with weak priors
#> that keep estimates off the boundary of the parameter space (Galindo
#> Garre & Vermunt, 2006): 1 pseudo-observation for the class proportions
#> and 1 pseudo-observation for the item-response probabilities. The
#> search used 20 random starts; the best log-likelihood was reached by 20
#> of them. The algorithm converged in 58 iterations.
#> 
#> Models with 1 to 4 classes were compared, including the one-class
#> model. The comparison reported the log-likelihood, AIC, BIC,
#> sample-size-adjusted BIC, CAIC, AIC3, AWE and relative entropy for
#> each. BIC was lowest for K = 4, the model reported here. Entropy is a
#> classification diagnostic, not a criterion for the number of classes
#> (Masyn, 2013; Sorgente et al., 2025).
#> 
#> The model's log-likelihood was -2531.14 with 67 free parameters (AIC =
#> 5196.28, BIC = 5463.71, sample-size-adjusted BIC = 5251.11).
#> 
#> Relative entropy was 0.907. The smallest class held 14.4% of the sample
#> (about 58 cases). The average posterior probability of the assigned
#> class ranged from 0.914 to 1.000 across classes. These describe how
#> cleanly cases are classified; they are not criteria for the number of
#> classes (Masyn, 2013).
#> 
#> Models were estimated with the R package mixtureEM (version 0.6.0) in R
#> 4.6.1. Random starts were drawn with random_state = 1. The models are
#> refitted by: compare_mixtures(X = items, k_range = 1:4, measurement =
#> "binary", n_init = 20, random_state = 1)
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
#> Sorgente, A., Caliciuri, R., Robba, M., Lanz, M., & Zumbo, B. D.
#>   (2025). A systematic review of latent class analysis in psychology:
#>   Examining the gap between guidelines and research practice. Behavior
#>   Research Methods, 57, 301.
```

Each paragraph is one element of a character vector, so you can paste
them into a manuscript one at a time, or edit them first. The references
are kept in `attr(, "references")` and printed once at the end, which
works for both author-year and numbered citation styles.

Two things the report tells you are worth acting on before you write:

- **Whether the best log-likelihood was replicated.** If the paragraph
  says it was reached by one start only, refit with more starts
  (`n_init = 100`) and check that nothing better appears
  ([`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md)).
- **Whether the sample you analysed is the one you collected.** Cases
  with no observed indicator are removed before estimation, and the
  first paragraph says how many.

To describe a different K than the BIC’s choice, pass `k`:
`report_methods(cmp, k = 3)`. The enumeration paragraph then says that
the reported model is not the lowest-BIC one.

## A trajectory model and the GRoLTS checklist

[`grolts_report()`](https://pdvalencia.github.io/mixtureEM/reference/grolts_report.md)
takes a growth model from
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
or a
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
comparison of them, and sorts the 21 GRoLTS items into three groups:
those the fit answers, those you must supply, and those that do not
apply to the model. Here is a small simulated example with two
trajectory classes over four occasions:

``` r

set.seed(5)
y <- matrix(rnorm(800, sd = 1), 200) +
  outer(rep(c(0, 1.5), each = 100), 0:3, function(a, t) a + (a / 2) * t)
traj <- compare_longitudinal(y, k_range = 1:3, model = "lcga", times = 4,
                             family = "gaussian", n_init = 5,
                             random_state = 1, verbose = FALSE)
grolts_report(traj)
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
#>     The search used 5 random starts; the best log-likelihood was
#>     reached by 5 of them. The algorithm converged in 9 iterations.
#> 
#> 10. Are the model comparison (and selection) tools described from a
#>     statistical perspective?
#>     Models with 1 to 3 classes were compared, including the one-class
#>     model. The comparison reported the log-likelihood, AIC, BIC,
#>     sample-size-adjusted BIC, CAIC, AIC3, AWE and relative entropy for
#>     each. BIC was lowest for K = 2, the model reported here. Entropy is
#>     a classification diagnostic, not a criterion for the number of
#>     classes (Masyn, 2013; Sorgente et al., 2025).
#> 
#> 11. Are the total number of fitted models reported, including a
#>     one-class solution?
#>     3 models were fitted, K = 1, 2, 3; the one-class model was
#>     included.
#> 
#> 12. Are the number of cases per class reported for each model?
#>     K = 1: 100.0%; K = 2: 50.6%, 49.4%; K = 3: 50.4%, 34.3%, 15.2%
#> 
#> 13. If classification of cases in a trajectory is the goal, is entropy
#>     reported?
#>     Relative entropy 0.980.
#> 
#> 15. Are characteristics of the final class solution numerically
#>     described?
#>     Growth-factor means by class are in summary(fit); class shares
#>     50.6%, 49.4%.
#> 
#> 16. Are the syntax files available?
#>     Models were estimated with the R package mixtureEM (version 0.6.0)
#>     in R 4.6.1. Random starts were drawn with random_state = 1. The
#>     models are refitted by: compare_longitudinal(indicators = y,
#>     k_range = 1:3, model = "lcga", times = 4, verbose = FALSE, family =
#>     "gaussian", n_init = 5, random_state = 1)
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
```

The “you must supply” group is the useful part. It lists what no fit can
know: how much measurement times varied within a wave, why the data are
missing, which other model shapes you tried and why you kept this one,
and the plots a reader needs to see.

The function prints no score such as “15 of 21”. A count treats the
items as interchangeable, and they are not; which gaps matter depends on
what the paper claims.

## References

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

Sorgente, A., Caliciuri, R., Robba, M., Lanz, M., & Zumbo, B. D. (2025).
A systematic review of latent class analysis in psychology: Examining
the gap between guidelines and research practice. *Behavior Research
Methods*, *57*, 301.

van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., &
Vermunt, J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on
latent trajectory studies. *Structural Equation Modeling*, *24*(3),
451-467. <https://doi.org/10.1080/10705511.2016.1247646>
