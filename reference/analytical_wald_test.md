# Analytical Wald Test for a Single Covariate

Performs an omnibus Wald chi-squared test for the effect of a single
covariate on latent class membership. The null hypothesis is that the
covariate has no effect on class membership across all non-reference
classes simultaneously; the degrees of freedom are (classes - 1) x
(covariate columns), the conventional df for this test.

The test uses whichever analytical variance the model was fitted with
(see
[`covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md))
and names it in the returned `Method` column. For small samples or
poorly conditioned Hessians, consider
[`wald_omnibus_test`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
instead, which uses bootstrapped standard errors from
[`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md).

A caveat that applies to any Wald test of a strong effect: under the
Hauck-Donner phenomenon a coefficient large enough to nearly separate a
class drives the Wald statistic back towards zero. A non-significant
omnibus test standing beside large coefficients is therefore a warning,
not a green light; a likelihood-ratio or bootstrap test settles it.

## Usage

``` r
analytical_wald_test(model, term_name, ref_class = 1)
```

## Arguments

- model:

  A fitted `mixture_model` object with a covariate structural model
  (fitted with `structural = "covariate"`).

- term_name:

  Character string. The name of the covariate to test, as supplied to
  the model. For a factor this is the variable, not one of its dummy
  columns: naming a three-level `Marital` tests both of its contrasts
  jointly, which is the point of the omnibus test. Naming a single
  column (`"Marital.Single"`) tests that column alone. Models fitted
  before the term grouping was recorded fall back to substring matching
  against the column names, with a warning when that is ambiguous.

- ref_class:

  Integer. The reference latent class. Default is `1`.

## Value

A single-row data frame with columns `Covariate`, `Wald_Chi2`, `df`,
`p_value`, and `Method`, the name of the variance estimator behind the
statistic.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
Z <- matrix(rnorm(100), nrow = 100)
colnames(Z) <- "age"
fit <- fit_mixture(X, Y = Z, n_components = 3, measurement = "binary",
                   structural = "covariate",
                   n_steps = 3, correction = "ML", n_init = 5)
analytical_wald_test(fit, term_name = "age")
} # }
```
