# Bootstrap Wald Omnibus Test for a Covariate

Performs a multiparameter Wald chi-squared test for the effect of a
single covariate on latent class membership, using bootstrapped standard
errors from
[`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md).
Tests the joint null hypothesis that the covariate has no effect across
all non-reference classes simultaneously.

Compared to
[`analytical_wald_test`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md),
this test rests on no asymptotic variance approximation at all, and
carries the step-1 variability the analytical estimators approximate to
first order (see
[`covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md)),
so it is worth the wait in small samples or when classes overlap
substantially.

## Usage

``` r
wald_omnibus_test(boot_results, term_name, assume_independence = TRUE)
```

## Arguments

- boot_results:

  A list returned by
  [`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md).

- term_name:

  Character string. The name of the covariate to test. Matched as a
  substring against the column names of `boot_results$orig_betas` (i.e.,
  partial matches work).

- assume_independence:

  Logical. If `TRUE` (default), the Wald statistic is computed using
  only the diagonal of the bootstrap covariance matrix, treating
  parameters as independent. If `FALSE`, the full covariance matrix is
  used via the Moore-Penrose pseudoinverse, which captures correlations
  between parameters but requires more replications for stable
  estimates.

## Value

A single-row data frame with columns `Covariate`, `Wald_Chi2`, `df`, and
`p_value`.

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
boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
wald_omnibus_test(boot, term_name = "age")

# Using the full covariance matrix (requires more bootstrap reps)
wald_omnibus_test(boot, term_name = "age", assume_independence = FALSE)
} # }
```
