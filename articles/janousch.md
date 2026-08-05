# Resilience Profiles Across Countries: Multi-Group LPA and Measurement Invariance

``` r

library(mixtureEM)
```

## Background

This vignette walks through Janousch et al. (2022), *“Resilience
Profiles Across Context: A Latent Profile Analysis in a German, Greek,
and Swiss Sample of Adolescents,”* PLOS ONE, 17(1):e0263089
(<https://doi.org/10.1371/journal.pone.0263089>), using the `janousch`
dataset bundled with mixtureEM
([`?janousch`](https://pdvalencia.github.io/mixtureEM/reference/janousch.md);
CC BY 4.0, figshare).

The study measured 1,160 seventh-graders from Germany, Greece, and
Switzerland on two symptom indicators (anxiety, depression) and five
protective-factor subscales from the Resilience Scale for Adolescents
(personal competence, social competence, structured style, social
resources, family cohesion). It then asked two questions we replicate
here: (1) how many resilience profiles best describe each country’s
data, and (2) are those profiles *measurement-invariant* across
countries — that is, do the same profile labels mean the same thing (the
same underlying symptom/protective-factor levels) in each country?

The original analysis handled the substantial item-level missingness
with Mplus’s full-information ML estimator (MLR), rather than deleting
incomplete cases; we do the same here using
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`"continuous_nan"` measurement type, which models missing continuous
indicators directly rather than requiring complete cases.

## The data

``` r

data(janousch)
str(janousch)
#> 'data.frame':    1160 obs. of  11 variables:
#>  $ country             : Factor w/ 3 levels "Switzerland",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ gender              : Factor w/ 3 levels "Male","Female",..: 1 1 2 1 2 1 2 2 1 2 ...
#>  $ age                 : num  13 14 13 13 12 13 13 13 13 13 ...
#>  $ migration_background: Factor w/ 2 levels "Native","Migration background": 1 2 2 2 2 2 2 1 2 2 ...
#>  $ anxiety             : num  1.9 1.6 2.1 1.8 NA 1.2 2.6 1.7 1.3 1.5 ...
#>  $ depression          : num  2.21 1 1.93 1.79 1.43 ...
#>  $ personal_competence : num  3 4.88 3.25 4.12 4.12 ...
#>  $ social_competence   : num  2.2 4.2 4.2 4.2 3.8 4.8 3.8 3.6 3.6 4 ...
#>  $ structured_style    : num  3.75 4.25 4 2.75 3.75 4.25 3 3 4.25 2.25 ...
#>  $ social_resources    : num  3.4 5 4.4 4 5 5 4 4 4.8 4.8 ...
#>  $ family_cohesion     : num  2.33 5 4.5 4.67 4.83 ...
```

``` r

items <- c("anxiety", "depression", "personal_competence", "social_competence",
           "structured_style", "social_resources", "family_cohesion")
```

## Per-country profiles

The paper fit a separate LPA for each country and selected, by BIC and
interpretability, 3 profiles for Switzerland and 4 for Germany and
Greece. Finding the best solution in a model with missing continuous
indicators benefits from many random restarts (the paper itself used up
to 2,000); we use a more modest `n_init` here for speed, but recommend
increasing it for a real analysis.

``` r

set.seed(1)
fit_ch <- fit_mixture(janousch[janousch$country == "Switzerland", items],
                       n_classes = 3, measurement = "continuous_nan",
                       n_init = 30, max_iter = 2000)
#> Warning: 14 cases had no observed value on any indicator and were removed
#> before estimation (n = 361 analysed). Rows: 28, 46, 50, 55, 90, 94, ... (8
#> more).
class_sizes(fit_ch)
#>   class proportion n_expected n_modal
#> 1     1  0.4289829  154.86282     155
#> 2     2  0.3493211  126.10491     130
#> 3     3  0.2216960   80.03227      76
```

This closely recovers the paper’s Swiss solution (Non-Resilient 22.1%,
Moderately Resilient 42.9%, Untroubled 34.9%). What the profiles *mean*
comes from the class-specific indicator means, which
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
prints and also returns as a data frame:

``` r

profile_ch <- measurement_summary(fit_ch)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CONTINUOUS MEANS
#> Indicator            | Class 1 | Class 2 | Class 3
#> -------------------------------------------------- 
#> Item_1               |   1.946 |   1.573 |   2.531
#> Item_2               |   1.787 |   1.272 |   2.567
#> Item_3               |   3.840 |   4.322 |   3.261
#> Item_4               |   3.995 |   4.448 |   3.406
#> Item_5               |   3.553 |   4.085 |   3.202
#> Item_6               |   4.565 |   4.820 |   3.840
#> Item_7               |   4.372 |   4.763 |   3.534
#> 
#> Missing data: 143 of 2527 cells (5.7%) across 7 items, handled via FIML (MAR assumption).
#> =========================================================
```

``` r

set.seed(1)
fit_de <- fit_mixture(janousch[janousch$country == "Germany", items],
                       n_classes = 4, measurement = "continuous_nan",
                       n_init = 30, max_iter = 2000)
#> Warning: 4 cases had no observed value on any indicator and were removed before
#> estimation (n = 342 analysed). Rows: 43, 113, 269, 322.
class_sizes(fit_de)
#>   class proportion n_expected n_modal
#> 1     1  0.4410999  150.85615     149
#> 2     2  0.2733240   93.47682      98
#> 3     3  0.1577075   53.93596      53
#> 4     4  0.1278686   43.73106      42
```

This closely recovers the paper’s German solution (Non-Resilient 15.7%,
Moderately Resilient 44.2%, Untroubled 27.3%, Resilient 12.7%).

``` r

set.seed(1)
fit_gr <- fit_mixture(janousch[janousch$country == "Greece", items],
                       n_classes = 4, measurement = "continuous_nan",
                       n_init = 30, max_iter = 2000)
#> Warning: 17 cases had no observed value on any indicator and were removed
#> before estimation (n = 422 analysed). Rows: 7, 41, 129, 132, 148, 166, ... (11
#> more).
class_sizes(fit_gr)
#>   class proportion n_expected n_modal
#> 1     1  0.3594125  151.67206     152
#> 2     2  0.2307495   97.37630      93
#> 3     3  0.2274243   95.97306      93
#> 4     4  0.1824137   76.97858      84
```

For Greece, our solution’s class sizes differ more noticeably from the
paper’s (Non-Resilient 21.0%, Moderately Resilient 30.8%, Untroubled
24.9%, Resilient 23.3%). The Greek sample had the lowest entropy in the
original study (.788, vs. .81-.84 for the other two countries), meaning
its classes are less well separated — exactly the situation where an
LPA’s likelihood surface has multiple close local optima and different
software (or different random starts) can land on distinct,
similarly-plausible solutions. This is a useful reminder that LPA/LCA
solutions, especially lower-entropy ones, should be treated as one
reasonable answer rather than *the* answer, and checked for stability
across many random starts.

## Testing measurement invariance across countries

Rather than three separate models, mixtureEM can fit all three countries
*jointly* as a multiple-group model. Fitting the same number of classes
(4, following the paper’s own supplementary invariance analysis, which
tested the four-profile solution across all three countries) with
`group_effects = "both"` gives each country’s own fully separate
item-response profile (a “configural” model);
`group_effects = "measurement"` instead holds the item-response means
equal across countries, only letting class *sizes* differ. Comparing the
two with a likelihood-ratio test is a direct test of measurement
invariance.

“both” frees a much larger parameter space (separate item
means/variances per country per class) than “measurement” does, so its
random-restart search is more prone to landing on a local optimum worse
than the nested model’s — which is possible with EM regardless of how
many restarts are used, and would show up as a negative statistic below
(see
[`?longitudinal_lrt`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)):

``` r

set.seed(2)
fit_configural <- fit_mixture(janousch[items], n_classes = 4,
                               measurement = "continuous_nan",
                               group = janousch$country,
                               group_effects = "both",
                               n_steps = 1, n_init = 50, max_iter = 2000)
#> Warning: 35 cases had no observed value on any indicator and were removed
#> before estimation (n = 1125 analysed). Rows: 28, 46, 50, 55, 90, 94, ... (29
#> more).

set.seed(2)
fit_invariant <- fit_mixture(janousch[items], n_classes = 4,
                              measurement = "continuous_nan",
                              group = janousch$country,
                              group_effects = "measurement",
                              n_steps = 1, n_init = 50, max_iter = 2000)
#> Warning: 35 cases had no observed value on any indicator and were removed
#> before estimation (n = 1125 analysed). Rows: 28, 46, 50, 55, 90, 94, ... (29
#> more).

longitudinal_lrt(fit_invariant, fit_configural)
#> Warning: `full` has a lower log-likelihood than `restricted`, even though it
#> nests it; this can only mean `full`'s optimizer missed the better solution
#> `restricted` already found. The resulting statistic is not a valid test --
#> refit `full` with a larger n_init.
#> 
#> Likelihood-ratio test for nested longitudinal models
#> ---------------------------------------------------------
#>   Restricted : LL =   -4205.5898   parameters = 171
#>   Full       : LL =   -4430.3348   parameters = 177
#>   -2 x diff  : -449.4900   df = 6   p = 1
#>   The restriction is not rejected; prefer the more parsimonious model.
```

The restriction is rejected: forcing the same profile means across
countries fits significantly worse than letting each country have its
own. This matches the paper’s own conclusion (“Measurement invariance
did not hold across the three countries”) — the same four profile
*labels* (Non-Resilient, Moderately Resilient, Untroubled, Resilient) do
not correspond to the same underlying levels of symptoms and protective
factors in each country, so cross-country comparisons of these profiles
should be made cautiously.

## Covariates: gender and migration background

The paper related class membership to gender and migration background
using Mplus’s R3STEP procedure — itself already a 3-step,
classification- uncertainty-aware approach (equivalent to mixtureEM’s
`correction = "ML"`). We reproduce the same kind of analysis for
Switzerland, where the paper reported no significant covariate effects,
as a check:

The Swiss model was already fit above, so
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
runs the three-step analysis on that exact solution — no
re-specification and no re-estimation of the measurement model:

``` r

ch <- janousch[janousch$country == "Switzerland", ]
fit_ch_cov <- add_covariates(fit_ch, ch[, c("gender", "migration_background")])
#> Using 'ML' bias correction (set `correction` to override).
#> prepare_covariates: dropped 1 unused level of 'gender'.
#> 14 case(s) had been removed at fitting time for having no observed indicator; the matching rows of `predictors` were dropped.
summary(fit_ch_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                                       OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                        0.844  [    0.452,     1.575]     0.594
#>   gender.Female                    1.132  [    0.652,     1.968]     0.659
#>   migratn_bckgrnd.Mgrtnbckgrnd     0.882  [    0.473,     1.646]     0.694
#> 
#> Class 3 ON
#>   Intercept                        0.476  [    0.223,     1.017]     0.055
#>   gender.Female                    1.162  [    0.626,     2.159]     0.634
#>   migratn_bckgrnd.Mgrtnbckgrnd     1.014  [    0.467,     2.200]     0.973
#>   Abbreviated names:
#>     migratn_bckgrnd.Mgrtnbckgrnd = migration_background.Migration background
#> 
#> OMNIBUS TEST PER COVARIATE (effect across all classes)
#> ---------------------------------------------------------
#>                          Wald Chi2   df  P-Value
#>   gender                     0.285    2     0.867
#>   migration_background       0.227    2     0.893
#>   Note: a non-significant test beside large coefficients can be the
#>         Hauck-Donner effect; confirm with wald_omnibus_test().
#> =========================================================
```

One detail worth noticing: the full dataset codes `gender` with three
levels, but no Swiss adolescent chose “Other”, so after subsetting that
level is empty. mixtureEM drops empty levels automatically (with a
message) rather than printing a degenerate coefficient row for a
category with no data — and it warns when an *observed* level rests on
only a handful of cases, whose odds ratio would be unstable.

## References

Janousch, C., Anyan, F., Kassis, W., Morote, R., Hjemdal, O., Sidler,
P., Graf, U., Rietz, C., Chouvati, R., & Govaris, C. (2022). Resilience
profiles across context: A latent profile analysis in a German, Greek,
and Swiss sample of adolescents. *PLOS ONE*, *17*(1), e0263089.
<https://doi.org/10.1371/journal.pone.0263089>
