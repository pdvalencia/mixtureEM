# LCA with Complex Survey Data

``` r

library(mixtureEM)
```

Large public-health and education surveys are almost never simple random
samples: cases carry *sampling weights* (they stand for different
numbers of people), and they arrive in *strata* and *clusters* (schools,
PSUs) whose members resemble each other. Ignoring the weights biases the
estimates; ignoring the design makes the standard errors too small.
mixtureEM handles both through the `weights`, `strata`, and `cluster`
arguments, available in every model.

## The data: 12 health-risk behaviors, YRBS 2005

The bundled `yrbs2005` dataset holds twelve dichotomous risk-behavior
items for 13,840 U.S. high-school students from the CDC’s 2005 Youth
Risk Behavior Survey, together with the survey’s weight, PSU, and
stratum variables
([`?yrbs2005`](https://pdvalencia.github.io/mixtureEM/reference/yrbs2005.md))
— Collins and Lanza’s (2010, ch. 2 and 4) own analysis sample, after
dropping cases missing on `grade` or on every item. We follow their
five-class solution.

``` r

data(yrbs2005)
items <- yrbs2005[, 3:14]
names(items)
#>  [1] "smoked_before_13"      "smoked_daily_30d"      "drove_drinking"       
#>  [4] "first_drink_before_13" "binge_drink_30d"       "marijuana_before_13"  
#>  [7] "cocaine_ever"          "glue_ever"             "meth_ever"            
#> [10] "ecstasy_ever"          "sex_before_13"         "sex_4plus_partners"
```

The items still contain some missing responses; mixtureEM switches
affected indicators to a full-information (FIML) estimator
automatically.

## A design-based five-class LCA

``` r

set.seed(1)
fit <- fit_mixture(items, n_classes = 5, measurement = "binary",
                   weights = yrbs2005$weight,
                   strata  = yrbs2005$stratum,
                   cluster = yrbs2005$psu,
                   n_init = 20, max_iter = 2000)
fit
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 5
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 502 iterations)
#> Missing Data       : 7186 / 166080 cells (4.3%) in 12 items — FIML (MAR assumption)
#> ---------------------------------------------------------
#>   Log-Likelihood : -47487.93
#>   Parameters     : 64
#>   AIC            : 95103.86
#>   BIC            : 95586.12
#>   SABIC          : 95382.73
#>   Rel. Entropy   : 0.8343
#>   Best solution  : found by 5 of 20 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 70.24%
#>   Class 2: 12.59%
#>   Class 3: 8.26%
#>   Class 4: 4.91%
#>   Class 5: 4.00%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```

The five classes reproduce the familiar YRBS typology — a large low-risk
class, substance-experimentation classes, and a small high-everything
class:

``` r

plot(fit, class_labels = c("Low risk", "Early experimenters",
                           "Alcohol-focused", "Hard drug users",
                           "High risk"),
     main = "Health-risk behavior classes (YRBS 2005)")
```

![](survey_lca_files/figure-html/plot-1.png)

These five classes are recognisably the YRBS typology, but they are not
identical to the solution Collins and Lanza report. They fitted the
items unweighted. The fit above uses the sampling weights, so it
estimates the class structure of the *population* of U.S. high-school
students rather than of the achieved sample, and the two are not the
same when the design oversamples some students and undersamples others.
Neither solution is a mistake; they answer slightly different questions,
and it is worth being explicit about which one a published table is
reporting.

``` r

params <- measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator             | Overall | Class 1 | Class 2 | Class 3 | Class 4 | Class 5
#> --------------------------------------------------------------------------------- 
#> smoked_before_13      |   0.159 |   0.039 |   0.614 |   0.132 |   0.206 |   0.866
#> smoked_daily_30d      |   0.134 |   0.026 |   0.301 |   0.252 |   0.558 |   0.734
#> drove_drinking        |   0.099 |   0.009 |   0.107 |   0.537 |   0.293 |   0.487
#> first_drink_before_13 |   0.255 |   0.138 |   0.693 |   0.238 |   0.205 |   0.894
#> binge_drink_30d       |   0.255 |   0.085 |   0.428 |   0.969 |   0.602 |   0.857
#> marijuana_before_13   |   0.086 |   0.004 |   0.371 |   0.038 |   0.058 |   0.774
#> cocaine_ever          |   0.076 |   0.003 |   0.042 |   0.056 |   0.644 |   0.843
#> glue_ever             |   0.124 |   0.061 |   0.185 |   0.148 |   0.416 |   0.632
#> meth_ever             |   0.061 |   0.003 |   0.023 |   0.016 |   0.564 |   0.677
#> ecstasy_ever          |   0.062 |   0.005 |   0.063 |   0.067 |   0.412 |   0.634
#> sex_before_13         |   0.062 |   0.020 |   0.243 |   0.011 |   0.035 |   0.364
#> sex_4plus_partners    |   0.141 |   0.054 |   0.303 |   0.282 |   0.380 |   0.579
#> 
#> Missing data: 7186 of 166080 cells (4.3%) across 12 items, handled via FIML (MAR assumption).
#> =========================================================
```

## Predictors under the survey design

The stored design travels with the fitted model, so the stepwise
covariate analysis is design-based too: the standard errors below are
linearized (sandwich) estimates that respect the strata and clusters.

``` r

fit_cov <- add_covariates(fit, yrbs2005[, "sex", drop = FALSE])
#> Using 'ML' bias correction (set `correction` to override).
results <- summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (survey-linearized step 3, hessian step 1)
#> ---------------------------------------------------------
#>                               OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                0.126  [    0.103,     0.155]    < .001
#>   sex.Male                 1.876  [    1.561,     2.255]    < .001
#> 
#> Class 3 ON
#>   Intercept                0.106  [    0.081,     0.140]    < .001
#>   sex.Male                 1.217  [    0.986,     1.502]     0.068
#> 
#> Class 4 ON
#>   Intercept                0.082  [    0.065,     0.104]    < .001
#>   sex.Male                 0.685  [    0.545,     0.861]     0.001
#> 
#> Class 5 ON
#>   Intercept                0.037  [    0.029,     0.047]    < .001
#>   sex.Male                 2.124  [    1.641,     2.748]    < .001
#> 
#> OMNIBUS TEST PER COVARIATE (effect across all classes)
#> ---------------------------------------------------------
#>                          Wald Chi2   df  P-Value
#>   sex                       75.168    4    < .001
#>   Note: a non-significant test beside large coefficients can be the
#>         Hauck-Donner effect; confirm with wald_omnibus_test().
#> =========================================================
```

Note the `Standard errors:` line naming the estimator: with a cluster
sample, a model that ignored the design would report noticeably narrower
intervals — narrower than the data warrant.

## What the weights do and do not change

Two practical points that trip up survey LCA:

- **Weight scale.** Sampling weights are rescaled internally so that
  they sum to the number of observed cases; only their *relative* sizes
  matter, and the effective sample size in BIC stays the number of cases
  actually collected. Frequency weights (`weight_type = "frequency"`,
  where one row stands for that many identical cases) are used as-is.
- **Class enumeration.** Information criteria remain comparable across
  class counts under a fixed design, so
  [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
  works unchanged — pass the same `weights`/`strata`/`cluster` to every
  candidate model.

## Should you weight at all?

Weighting is not free. Weights protect against bias from unequal
selection probabilities and from coverage and nonresponse error, but a
weighted estimator is less efficient than an unweighted one; if the
weights were not needed, all you have bought is a wider confidence
interval. Bollen et al. (2016) observe that researchers rarely test
whether the weights are necessary and are guided more by disciplinary
convention than by evidence.

The diagnostics they review fall into two families: *difference in
coefficients* tests, which compare the weighted and unweighted estimates
and ask whether the gap is larger than sampling error (a Hausman-type
comparison), and *weight association* tests, which ask whether the
weight itself still predicts the outcome once the model’s covariates are
in. The two families are closely related, and the review is candid that
their finite-sample properties are not fully settled.

mixtureEM does not implement such a test, and there is no accepted
version of one for a latent class model. What it does make cheap is the
informal version: fit the same model with and without `weights =` and
compare the class sizes and item probabilities. Doing that here, the
weighted and unweighted five-class solutions do not just shift class
sizes by a point or two — several classes reorder which items they are
elevated on (the unweighted fit’s smallest classes read very differently
from the weighted fit’s, for instance on `drove_drinking` and
`sex_before_13`). The weights are carrying real information about who is
and is not in the achieved sample, and dropping them here would bias the
solution rather than merely widen its standard errors.

## References

Bollen, K. A., Biemer, P. P., Karr, A. F., Tueller, S., & Berzofsky, M.
E. (2016). Are survey weights needed? A review of diagnostic tests in
regression analysis. *Annual Review of Statistics and Its Application*,
3, 375-392. <https://doi.org/10.1146/annurev-statistics-011516-012958>.

Collins, L. M., & Lanza, S. T. (2010). *Latent class and latent
transition analysis: With applications in the social, behavioral, and
health sciences*. Wiley.
