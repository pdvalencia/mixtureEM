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

Applied papers usually name a growth mixture model by what the classes
are allowed to differ in, and Ram and Grimm (2009, pp. 569-570) give the
three levels. *Means* — classes differ in mean change only — is
`psi = "equal"`, the default here. *Means + Covs* — classes also differ
in how much interindividual variation there is within them — is
`psi = "free"`; the residual variances are a separate switch,
`residual_equal = FALSE`, and should not be merged with it. *Means +
Covs + Pattern*, where classes differ in the shape of change as well,
would need class-specific time scores; `fit_gmm()` builds one design for
every class, so that level is not available.

Two cautions worth carrying into any such analysis. Classes will be
found whether or not there are groups to find: a skewed or otherwise
non-normal outcome can be fitted by extra classes that are not subgroups
of anyone (Jung & Wickrama, 2008, p. 305; Lee et al., 2023, p. 652, both
citing Bauer & Curran, 2003). And as Ram and Grimm put it (p. 574),
"groups and differences among groups will be found; but whether they
represent true processes that generated the data is unknown" — the
remedies they name are replication on new data and checking that class
membership relates to other measured variables in the ways theory
predicts.

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
  (intercept and slope), `2` for a quadratic, and so on. Each degree
  needs occasions to identify it — with three time points a linear
  pattern can be modelled, with four a quadratic as well, and with five
  a cubic (Berlin et al., 2014, p. 191). Asking for more than the
  occasions support is an error rather than a warning.

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

  It is not free, though, and the cost runs the other way. Simulation
  work reported by Lee et al. (2023, p. 651) finds that equality
  restrictions on the growth-factor variances and covariances "could
  result in the over-extraction of latent classes and biased parameter
  estimates" — the within-class heterogeneity the constraint refuses to
  let differ between classes has to go somewhere, and an extra class is
  where it goes. So `psi = "equal"` buys stability and bounds the
  likelihood, and it can buy an extra class that is an artefact of the
  constraint. Where the number of classes is itself the finding, fit
  both and report which was used.

  This pulls against the collapsed-variance warning, which fires more
  often under `psi = "free"`. That is not a contradiction: they are the
  two sides of one trade-off between a model flexible enough to be
  realistic and one constrained enough to be estimable.

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
#>   Parameters     : 12
#>   AIC            : 5221.01
#>   BIC            : 5268.90
#>   SABIC          : 5230.83
#>   Rel. Entropy   : 0.6348
#>   Best solution  : found by 2 of 4 starts that ran to convergence (of 10 requested)
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 54.71%
#>   Class 2: 45.29%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
# }
```
