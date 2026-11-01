# Constructor for structured multivariate-normal emissions

Sets up the emission state for a growth mixture model: each class
follows its own polynomial trajectory in the growth-factor means, and
cases vary about their class's trajectory through random effects with a
class-specific or class-common covariance.

## Usage

``` r
structured_normal_model(
  n_components,
  design,
  random_effects = "intercept_slope",
  psi = "equal",
  residual = "occasion",
  residual_equal = TRUE,
  growth_covariates = NULL,
  growth_covariates_equal = TRUE,
  psi_fixed = NULL,
  theta_pairs = NULL,
  ...
)
```

## Arguments

- n_components:

  Integer. The number of latent classes to estimate.

- design:

  Numeric matrix with one row per occasion and one column per growth
  coefficient, as built by the polynomial in the time scores.

- random_effects:

  Which growth factors vary within a class: `"none"`, `"intercept"`,
  `"intercept_slope"`, or `"all"`; or the columns of `design` that do,
  as integers.

- psi:

  Whether the growth-factor covariance is held `"equal"` across classes,
  has class-specific `"variances"` with the covariances held equal, or
  is estimated `"free"`ly in each.

- residual:

  Whether residual variances are `"occasion"`-specific or `"constant"`
  across occasions.

- residual_equal:

  Logical. Hold the residual variances equal across classes.

- growth_covariates:

  Optional numeric matrix of case-level covariates, one row per case and
  one column per covariate, regressing the growth factors on external
  variables. Must be row-aligned with the outcome matrix the emission is
  fitted to, and complete.

- growth_covariates_equal:

  Logical. Hold those regressions equal across classes.

- psi_fixed:

  Optional list of one matrix per class, the size of \\\Psi\\, holding
  the values fixed in that class and `NA` where the element is
  estimated. Built by
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  from its `fix` argument.

- theta_pairs:

  Optional two-column integer matrix, one row per pair of columns of the
  outcome whose residuals covary, each pair once. The residual
  covariance is then \\\Theta\\ with those cells free, shared by the
  classes or not as `residual_equal` says. Built by
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  from its `residual_covariances` argument.

- ...:

  Additional arguments, ignored.

## Value

A list object of class `c("structured_normal", "emission")`.
