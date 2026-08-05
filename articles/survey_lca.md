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
items for 13,917 U.S. high-school students from the CDC’s 2005 Youth
Risk Behavior Survey, together with the survey’s weight, PSU, and
stratum variables
([`?yrbs2005`](https://pdvalencia.github.io/mixtureEM/reference/yrbs2005.md)).
Collins and Lanza (2010, ch. 2 and 4) built their textbook LCA
illustration on these items; we follow their five-class solution.

``` r

data(yrbs2005)
items <- yrbs2005[, 3:14]
names(items)
#>  [1] "smoked_before_13"      "smoked_daily_30d"      "drove_drinking"       
#>  [4] "first_drink_before_13" "binge_drink_30d"       "marijuana_before_13"  
#>  [7] "cocaine_ever"          "glue_ever"             "meth_ever"            
#> [10] "ecstasy_ever"          "sex_before_13"         "sex_4plus_partners"
```

The items contain missing responses; mixtureEM switches affected
indicators to a full-information (FIML) estimator automatically, and a
handful of cases missing on *every* item are dropped with a warning.

## A design-based five-class LCA

``` r

set.seed(1)
fit <- suppressWarnings(
  fit_mixture(items, n_classes = 5, measurement = "binary",
              weights = yrbs2005$weight,
              strata  = yrbs2005$stratum,
              cluster = yrbs2005$psu,
              n_init = 20, max_iter = 2000)
)
fit
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 5
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 12 iterations)
#> Cases Removed      : 3 of 13917 with no observed indicator (n = 13914 analysed)
#> Missing Data       : 7284 / 166968 cells (4.4%) in 12 items — FIML (MAR assumption)
#> ---------------------------------------------------------
#>   Log-Likelihood : -47814.65
#>   Rel. Entropy   : 0.8340
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 70.12%
#>   Class 2: 12.68%
#>   Class 3: 8.25%
#>   Class 4: 4.92%
#>   Class 5: 4.03%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
class_sizes(fit)
#>   class proportion n_expected    n_modal
#> 1     1 0.70117440  9756.1407 10136.8160
#> 2     2 0.12683859  1764.8322  1592.7332
#> 3     3 0.08246433  1147.4087   997.4884
#> 4     4 0.04923695   685.0829   639.3167
#> 5     5 0.04028572   560.5355   547.6457
```

The five classes reproduce the familiar YRBS typology — a large low-risk
class, substance-experimentation classes, and a small high-everything
class:

``` r

plot(fit, main = "Health-risk behavior classes (YRBS 2005)")
```

![](survey_lca_files/figure-html/plot-1.png)

``` r

params <- measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator             | Class 1 | Class 2 | Class 3 | Class 4 | Class 5
#> ----------------------------------------------------------------------- 
#> smoked_before_13      |   0.039 |   0.615 |   0.130 |   0.207 |   0.871
#> smoked_daily_30d      |   0.026 |   0.302 |   0.251 |   0.561 |   0.737
#> drove_drinking        |   0.009 |   0.109 |   0.533 |   0.295 |   0.493
#> first_drink_before_13 |   0.138 |   0.695 |   0.236 |   0.207 |   0.891
#> binge_drink_30d       |   0.085 |   0.429 |   0.968 |   0.603 |   0.862
#> marijuana_before_13   |   0.005 |   0.370 |   0.037 |   0.059 |   0.780
#> cocaine_ever          |   0.003 |   0.042 |   0.056 |   0.645 |   0.846
#> glue_ever             |   0.060 |   0.186 |   0.147 |   0.417 |   0.636
#> meth_ever             |   0.003 |   0.024 |   0.016 |   0.565 |   0.680
#> ecstasy_ever          |   0.005 |   0.064 |   0.066 |   0.416 |   0.637
#> sex_before_13         |   0.020 |   0.245 |   0.011 |   0.034 |   0.376
#> sex_4plus_partners    |   0.054 |   0.306 |   0.282 |   0.382 |   0.588
#> 
#> Missing data: 7284 of 166968 cells (4.4%) across 12 items, handled via FIML (MAR assumption).
#> =========================================================
```

## Predictors under the survey design

The stored design travels with the fitted model, so the stepwise
covariate analysis is design-based too: the standard errors below are
linearized (sandwich) estimates that respect the strata and clusters.

``` r

fit_cov <- add_covariates(fit, yrbs2005[, c("grade", "sex")])
#> Using 'ML' bias correction (set `correction` to override).
#> 3 case(s) had been removed at fitting time for having no observed indicator; the matching rows of `predictors` were dropped.
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
#>   Intercept                0.180  [    0.146,     0.222]    < .001
#>   grade.10                 0.742  [    0.622,     0.885]    < .001
#>   grade.11                 0.574  [    0.448,     0.734]    < .001
#>   grade.12                 0.475  [    0.360,     0.626]    < .001
#>   sex.Male                 1.795  [    1.508,     2.136]    < .001
#> 
#> Class 3 ON
#>   Intercept                0.005  [    0.000,     0.124]     0.001
#>   grade.10                13.245  [    0.547,   320.836]     0.112
#>   grade.11                31.785  [    1.336,   756.475]     0.032
#>   grade.12                53.759  [    2.121,  1362.381]     0.016
#>   sex.Male                 1.310  [    1.061,     1.618]     0.012
#> 
#> Class 4 ON
#>   Intercept                0.042  [    0.026,     0.067]    < .001
#>   grade.10                 1.669  [    0.868,     3.209]     0.125
#>   grade.11                 2.900  [    1.878,     4.478]    < .001
#>   grade.12                 2.728  [    1.749,     4.256]    < .001
#>   sex.Male                 0.708  [    0.562,     0.893]     0.004
#> 
#> Class 5 ON
#>   Intercept                0.041  [    0.030,     0.056]    < .001
#>   grade.10                 0.928  [    0.637,     1.351]     0.696
#>   grade.11                 0.771  [    0.519,     1.144]     0.197
#>   grade.12                 0.912  [    0.651,     1.276]     0.590
#>   sex.Male                 2.125  [    1.640,     2.753]    < .001
#> 
#> OMNIBUS TEST PER COVARIATE (effect across all classes)
#> ---------------------------------------------------------
#>                          Wald Chi2   df  P-Value
#>   grade                    325.934   12    < .001
#>   sex                       70.689    4    < .001
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

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent class and latent
transition analysis: With applications in the social, behavioral, and
health sciences*. Wiley.
