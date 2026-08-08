# Growth Mixture Modeling

Fits class-specific growth curves to a continuous outcome measured
repeatedly, allowing cases to vary about their own class's trajectory
through random growth factors. This is the growth mixture model of
Muthen & Shedden (1999) and, with `random_effects = "none"`, latent
class growth analysis for a continuous outcome.

The model is \$\$y_i \mid c_i = k \sim N(\Lambda (\alpha_k + \Gamma_k
x_i),\\ \Lambda_r \Psi_k \Lambda_r' + \Theta_k),\$\$ where \\\Lambda\\
is the polynomial design in time, \\\alpha_k\\ the growth-factor means
of class \\k\\, \\\Psi_k\\ the covariance of the growth factors that are
allowed to vary within a class, and \\\Theta_k\\ the diagonal matrix of
residual variances. Because the outcome is continuous the random effects
integrate out in closed form, so no numerical integration is involved
and the estimator is the ordinary EM algorithm the rest of the package
uses.

\\\Gamma_k\\ is present only when `growth_predictors` are supplied, and
is the within-class growth-factor regression: covariates shift a case's
own growth factors *within* its class, answering "who, in this
trajectory group, starts higher or grows faster?". It is a different
question from `predictors`, which is `c ON x` and asks who is *in* the
group; the two can be asked together, and are then estimated in one pass
with `n_steps = 1`. With covariates present \\\alpha_k\\ is the
growth-factor *intercept* rather than its mean, and the reported
trajectory is evaluated at the sample mean of the covariates.

Contrast
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
which fixes the growth-factor variances at zero, so every case in a
class sits on that class's curve up to occasion-level noise. LCGA is the
more parsimonious model and, in criminology especially, often the
preferred one; GMM is the more realistic and is what developmental
psychology and epidemiology usually publish. The two are nested, so the
choice can be made on BIC.

## Usage

``` r
fit_gmm(
  indicator,
  n_classes = 2,
  times = NULL,
  degree = 1,
  random_effects = c("intercept_slope", "intercept", "none", "all"),
  psi = c("equal", "free"),
  residual = c("occasion", "constant"),
  residual_equal = TRUE,
  time_scores = NULL,
  layout = c("time_major", "item_major"),
  id = NULL,
  time = NULL,
  item = NULL,
  time_labels = NULL,
  predictors = NULL,
  growth_predictors = NULL,
  growth_predictors_equal = TRUE,
  ...
)
```

## Arguments

- indicator:

  The repeated outcome. Either a wide matrix or data frame with one
  column per occasion, a three-dimensional array with dimensions n by 1
  by times, or a long data frame together with `id` and `time`. Exactly
  one variable may be modelled.

- n_classes:

  Integer. Number of trajectory classes.

- times:

  Integer. Number of occasions. Required for wide input; inferred
  otherwise.

- degree:

  Degree of the polynomial in time: `1` for a linear trajectory
  (intercept and slope), `2` for a quadratic, and so on.

- random_effects:

  Which growth factors vary within a class:

  - `"intercept_slope"` (default) — both the intercept and the linear
    slope, the standard growth mixture model. Cases in a class differ
    both in where they start and in how fast they change.

  - `"intercept"` — the intercept only. Cases in a class start at
    different levels but change in parallel.

  - `"none"` — no random effects, which is latent class growth analysis;
    the same model as `fit_lcga(family = "gaussian")`, but with the
    residual variances free over occasions by default.

  - `"all"` — every growth factor, including a quadratic or higher term.

- psi:

  Whether the growth-factor covariance \\\Psi\\ is held `"equal"` across
  classes (the default) or estimated `"free"`ly in each. Equality is the
  better default in practice as well as by convention: it is what BIC
  typically selects, and a class-specific \\\Psi\\ is where growth
  mixture models most often produce a negative variance.

- residual:

  Whether the residual variances are `"occasion"`-specific (the default)
  or `"constant"` across occasions.

- residual_equal:

  Logical. Hold the residual variances equal across classes (the
  conventional default).

- time_scores:

  Numeric values of time used in the polynomial, one per occasion.
  Defaults to `0, 1, ..., times - 1`, which makes the intercept growth
  factor the level at the first occasion. Supply the actual measurement
  times when the occasions are unequally spaced.

- layout:

  For wide input with a three-dimensional array or several columns per
  occasion, whether columns run `"time_major"` or `"item_major"`.

- id, time:

  For long input, the case and occasion identifiers, given either as
  column names or as vectors.

- item:

  For long input, the column holding the outcome.

- time_labels:

  Optional display labels for the occasions.

- predictors:

  Optional predictors of class membership (`c ON x`), passed to
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)'s
  three-step machinery.

- growth_predictors:

  Optional covariates regressing the growth factors themselves
  (`i s ON x`): a vector, matrix or data frame with one row per case,
  factors dummy-coded. These enter the measurement model rather than the
  structural one, so they are estimated jointly with the trajectories in
  a single pass regardless of `n_steps`. They must be complete: a case
  missing a covariate has no model-implied mean, and unlike a missing
  occasion there is nothing here to integrate it out, so remove or
  impute those cases first.

- growth_predictors_equal:

  Logical. Hold the growth-factor regressions equal across classes (the
  conventional default). Set `FALSE` to let each class have its own,
  which asks whether a covariate matters differently in different
  trajectory groups.

- ...:

  Further arguments passed to
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  such as `outcome` and its companions, `n_init`, `random_state`,
  `weights`, `strata` or `cluster`.

## Value

An object of class `c("gmm", "mixture_model")`. In addition to the usual
fields it carries `$growth`, holding the design matrix, time scores, the
per-class growth-factor means, the fitted class trajectories, the
growth-factor covariance of each class, the residual variances and, when
`growth_predictors` were given, `$growth$coefficients`: one
growth-factor-by-covariate matrix per class.

## References

Muthen, B., & Shedden, K. (1999). Finite mixture modeling with mixture
outcomes using the EM algorithm. *Biometrics*, *55*(2), 463-469.

## See also

[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
for the no-random-effects version and
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
for trajectory classes with no growth curve at all.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 400
cls <- rbinom(n, 1, 0.5)
icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
fit <- fit_gmm(y, times = 4, n_classes = 2, n_init = 10, random_state = 1)
fit
#> 
#> =========================================================
#>              GROWTH MIXTURE MODEL
#> =========================================================
#> Occasions          : 4 (time scores 0, 1, 2, 3)
#> Trajectory         : linear
#> Random effects     : intercept, linear
#> 
#> GROWTH FACTOR MEANS
#>         intercept linear
#> Class 1     1.026  0.179
#> Class 2     2.950  1.003
#> 
#> GROWTH FACTOR (CO)VARIANCE (held equal across classes)
#>           intercept linear
#> intercept     1.204  0.005
#> linear        0.005  0.189
#> 
#> RESIDUAL VARIANCE (held equal across classes)
#>             T1    T2    T3    T4
#> Variance 0.509 0.537 0.572 0.336
#> 
#> FITTED TRAJECTORY (mean)
#>            T1    T2    T3    T4
#> Class 1 1.026 1.205 1.384 1.563
#> Class 2 2.950 3.952 4.955 5.958
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 125 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -2598.50
#>   Rel. Entropy   : 0.6348
#>   Best solution  : found by 2 of 4 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 54.71%
#>   Class 2: 45.29%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
# }
```
