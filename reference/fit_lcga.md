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

**What the zero-variance assumption costs when you enumerate classes.**
In an LCGA the information criteria can keep improving as classes are
added, because the within-class variation the model refuses to estimate
has to be absorbed somewhere, and more classes is where it goes. Berlin
et al. (2014, p. 196) hit this in their own analysis — "increasing the
number of latent classes resulted in increasingly better (i.e., smaller)
AICs, BICs, and SSA-BICs, without any detriment to entropy" — and read
it as a symptom rather than a result: "this is probably a poorly
specified model … The assumption of zero variance within classes (as
modeled here) is not likely tenable and might account for each
successive model … seeming to improve the fit." So before reading a
monotone BIC as evidence for many classes, fit a
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
at the same K and compare.

That is not an argument against the model. Where within-class
homogeneity genuinely holds, LCGAs "often allow for more straightforward
interpretations" (Ram & Grimm, 2009, p. 574), and the constraint "may be
particularly helpful when working with smaller sample sizes or when more
complex models fail due to nonconvergence, out of range estimates, or
other statistical problems, or as an initial modeling step prior to
specifying a GMM model" (Berlin et al., 2014, p. 191) — which is also
how Jung and Wickrama (2008, p. 304) recommend starting.

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
  `degree` cannot exceed `times - 2`. In practice that means three time
  points support a linear pattern, four a quadratic as well, and five a
  cubic (Berlin et al., 2014, p. 191).

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

## References

Nagin, D. S. (1999). Analyzing developmental trajectories: A
semiparametric, group-based approach. *Psychological Methods*, *4*(2),
139-157.
[doi:10.1037/1082-989X.4.2.139](https://doi.org/10.1037/1082-989X.4.2.139)

Nagin, D. S., & Odgers, C. L. (2010). Group-based trajectory modeling in
clinical research. *Annual Review of Clinical Psychology*, *6*, 109-138.
[doi:10.1146/annurev.clinpsy.121208.131413](https://doi.org/10.1146/annurev.clinpsy.121208.131413)

Berlin, K. S., Parra, G. R., & Williams, N. A. (2014). An introduction
to latent variable mixture modeling (part 2): Longitudinal latent class
growth analysis and growth mixture models. *Journal of Pediatric
Psychology*, *39*(2), 188-203.
[doi:10.1093/jpepsy/jst085](https://doi.org/10.1093/jpepsy/jst085)

Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class
growth analysis and growth mixture modeling. *Social and Personality
Psychology Compass*, *2*(1), 302-317.
[doi:10.1111/j.1751-9004.2007.00054.x](https://doi.org/10.1111/j.1751-9004.2007.00054.x)

Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction
to growth mixture models (GMM). In *International Encyclopedia of
Education* (4th ed., Vol. 14, pp. 646-655). Elsevier.
[doi:10.1016/B978-0-12-818630-5.10076-4](https://doi.org/10.1016/B978-0-12-818630-5.10076-4)

Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
identifying differences in longitudinal change among unobserved groups.
*International Journal of Behavioral Development*, *33*(6), 565-576.
[doi:10.1177/0165025409343765](https://doi.org/10.1177/0165025409343765)

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
#>   Best solution  : found by 5 of 5 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 58.21%
#>   Class 2: 41.79%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
# }
```
