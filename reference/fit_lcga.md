# Latent Class Growth Analysis

Fits class-specific growth curves to a single outcome measured
repeatedly: each latent class follows its own polynomial trajectory over
time, with no within-class random effects. This is Nagin's group-based
trajectory modelling (Nagin, 2005). Because there are no random effects,
the occasions are conditionally independent given class, so everything
the package already offers - model selection, the bootstrap
likelihood-ratio test, predictors of class membership, distal outcomes,
survey designs and FIML for missing occasions - applies unchanged.

The trajectory is modelled on the link scale, \$\$g(E\[y\_{it} \mid c_i
= k\]) = \beta\_{k0} + \beta\_{k1} t + \beta\_{k2} t^2 + \dots,\$\$ with
the logit link for a binary outcome. Contrast
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md),
which also produces trajectory classes but leaves each occasion's
response parameter free rather than constraining it to a smooth curve;
LCGA is the more parsimonious model and the one that extrapolates.

## Usage

``` r
fit_lcga(
  indicator,
  n_classes = 2,
  times = NULL,
  degree = 1,
  family = c("binomial", "gaussian", "poisson"),
  time_scores = NULL,
  layout = c("time_major", "item_major"),
  id = NULL,
  time = NULL,
  item = NULL,
  time_labels = NULL,
  predictors = NULL,
  ...
)
```

## Arguments

- indicator:

  The repeated outcome. Either a wide matrix or data frame with one
  column per occasion, a three-dimensional array with dimensions n by 1
  by times, or a long data frame together with `id` and `time`. Exactly
  one variable may be modelled; for several at once see
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).
  Named for consistency with
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)'s
  `indicators`, and so that `outcome=` stays free for a distal outcome
  caused by the trajectory class.

- n_classes:

  Integer. Number of trajectory classes.

- times:

  Integer. Number of occasions. Required for wide input; inferred
  otherwise.

- degree:

  Degree of the polynomial in time: `1` for a linear trajectory
  (intercept and slope), `2` for a quadratic, and so on. Must leave at
  least one occasion beyond the coefficients being estimated, so
  `degree` cannot exceed `times - 2`.

- family:

  Distribution of the outcome given class, which also sets the link the
  trajectory is linear on:

  - `"binomial"` — a binary outcome in `{0, 1}`, logit link. The
    trajectory is a curve in the probability of the event.

  - `"gaussian"` — a continuous outcome, identity link. The trajectory
    is a curve in the mean, and each class also carries a residual
    variance, constant across occasions. This is LCGA in the sense of
    Nagin's censored-normal model without the censoring, and the
    no-random-effects special case of a growth mixture model.

  - `"poisson"` — a count outcome, log link. The trajectory is a curve
    in the event rate. Counts must be whole and non-negative;
    over-dispersed or zero-inflated counts are not yet modelled as such.

- time_scores:

  Numeric values of time used in the polynomial, one per occasion.
  Defaults to `0, 1, ..., times - 1`, which makes the intercept the
  fitted value at the first occasion. Supply the actual measurement
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

  Optional predictors of class membership, passed to
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)'s
  three-step machinery.

- ...:

  Further arguments passed to
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  such as `outcome` and its companions, `n_init`, `random_state`,
  `weights`, `strata` or `cluster`.

## Value

An object of class `c("lcga", "mixture_model")`. In addition to the
usual fields it carries `$growth`, holding the design matrix, time
scores, family, per-class coefficients, the fitted trajectories, and —
for `family = "gaussian"` — the per-class residual variance.

## See also

[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
for unconstrained trajectory classes and
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
for a model in which class membership itself changes.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 400
cls <- rbinom(n, 1, 0.5)
eta <- outer(ifelse(cls == 1, 0.0, -0.8), rep(1, 4)) +
  outer(ifelse(cls == 1, 1.1, -0.1), 0:3)
y <- matrix(rbinom(n * 4, 1, plogis(eta)), n, 4)
fit <- fit_lcga(y, times = 4, n_classes = 2, n_init = 5, random_state = 1)
fit
#> 
#> =========================================================
#>            LATENT CLASS GROWTH ANALYSIS
#> =========================================================
#> Occasions          : 4 (time scores 0, 1, 2, 3)
#> Trajectory         : linear, logit link
#> 
#> GROWTH COEFFICIENTS (link scale)
#>         intercept linear
#> Class 1    -0.732 -0.055
#> Class 2    -0.048  1.124
#> 
#> FITTED TRAJECTORY (probability)
#>            T1    T2    T3    T4
#> Class 1 0.325 0.313 0.301 0.290
#> Class 2 0.488 0.746 0.900 0.965
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 73 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -1032.42
#>   Rel. Entropy   : 0.6318
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 58.21%
#>   Class 2: 41.79%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
# }
```
