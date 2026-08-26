# Fit a Latent Mixture Model (LCA / LPA)

The core estimation function. Fits a latent class analysis (LCA) or
latent profile analysis (LPA) model using the EM algorithm. Optionally
fits a structural model (covariates or distal outcomes) using 1-, 2-, or
3-step estimation with optional bias correction.

## Usage

``` r
fit_mixture_internal(
  X,
  Y = NULL,
  n_components = 2,
  measurement = "binary",
  structural = NULL,
  n_steps = 1,
  correction = "none",
  assignment = c("proportional", "modal"),
  n_init = 20,
  max_iter = 1000,
  random_state = NULL,
  order_by_size = TRUE,
  weights = NULL,
  weight_type = c("sampling", "frequency"),
  strata = NULL,
  cluster = NULL,
  refine = TRUE,
  bayes_constants = NULL,
  warm_start = NULL,
  se = c("corrected", "robust", "hessian"),
  n_cores = 1L,
  ...
)
```

## Arguments

- X:

  A numeric matrix or data frame of indicator variables for the
  measurement model. Rows are observations; columns are items or
  variables.

- Y:

  Optional numeric matrix or data frame of outcome or covariate
  variables for the structural model. Must be provided when `structural`
  is not `NULL`.

- n_components:

  Positive integer. Number of latent classes (or profiles) to estimate.
  Default is `2`.

- measurement:

  Character string or named list specifying the measurement model type.
  Accepted strings: `"binary"` / `"bernoulli"`, `"categorical"` /
  `"multinoulli"`, `"continuous"` / `"gaussian_diag"`, `"gaussian"` /
  `"gaussian_unit"`, `"count"` / `"poisson"`. Missing values are handled
  automatically: any indicator column containing `NA` is estimated with
  a full-information (FIML) variant that masks the missing cells under a
  missing-at-random assumption, while complete columns use the faster
  complete-data estimator. A single specification (e.g. `"binary"`)
  therefore covers both complete and incomplete data, and the fitted
  object reports any missingness it found. Cases missing on *every*
  indicator are the exception: they carry no information about class
  membership, so they are deleted before estimation and reported by
  [`print()`](https://rdrr.io/r/base/print.html) and in
  `$missing_data$n_empty_rows`. The explicit `"*_nan"` forms (e.g.
  `"binary_nan"`, `"continuous_nan"`) remain accepted as aliases that
  force the missing-data variant. Pass a named list to specify a mixed
  (nested) measurement model with different variable types; each block's
  missing-data handling is resolved from the columns it governs. Default
  is `"binary"`.

- structural:

  Character string specifying the structural model type. One of
  `"covariate"`, `"distal_regression"`, `"distal_pooled"`,
  `"distal_continuous"`, `"distal_continuous_regression"`. Requires `Y`.
  Default is `NULL` (measurement model only).

- n_steps:

  Integer. Estimation approach: `1` for simultaneous 1-step, `2` for
  2-step, or `3` for bias-corrected 3-step. Default is `1`.

- correction:

  Character. Bias correction for 3-step estimation. One of `"none"`,
  `"BCH"`, or `"ML"`. Ignored when `n_steps` is not `3`. Default is
  `"none"`.

- assignment:

  Character. How step-one class membership is carried into the
  step-three correction. `"proportional"` (default) uses the full
  posterior class probabilities; `"modal"` hardens each case to its most
  likely class first. Ignored when `n_steps` is not `3`.

- n_init:

  Positive integer. Number of random restarts. The solution with the
  highest log-likelihood is retained. Default is `20`.

- max_iter:

  Positive integer. Maximum EM iterations per restart. Default is
  `1000`.

- random_state:

  Optional integer seed for reproducibility. Default is `NULL`.

- order_by_size:

  Logical. If `TRUE` (default), classes are sorted from largest to
  smallest after fitting.

- weights:

  Optional numeric vector of length `nrow(X)` for survey or case
  weights. Default is `NULL` (equal weights of 1).

- weight_type:

  Either `"sampling"` (default; weights are rescaled to sum to the
  number of cases) or `"frequency"` (weights are counts of identical
  cases and set the effective sample size).

- strata:

  Optional vector of stratum identifiers for complex survey designs.

- cluster:

  Optional vector of cluster identifiers for complex survey designs.

- refine:

  Logical. If `TRUE` (default), applies L-BFGS refinement after EM
  convergence to optimize the penalized maximum likelihood.

- bayes_constants:

  Optional named list of prior strengths (`latent`, `categorical`,
  `poisson`, `variances`), each defaulting to `1`. See
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- warm_start:

  Optional function of `(model_state, X, Y)` returning a starting model
  state for EM, or `NULL` to skip that start. Used by the group-varying
  measurement search to seed each fit from the pooled solution. `NULL`
  (default) uses only the usual random initializations.

- se:

  Character. How standard errors for a covariate (class-prediction)
  structural model are computed when `n_steps` is `2` or `3`.
  `"corrected"` (default) is the first-order corrected estimator of Bakk
  et al. (2014): the step-3 sandwich plus the variance propagated from
  step 1. `"robust"` keeps only the sandwich. `"hessian"` inverts the
  step-3 observed information alone. See
  [`covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md)
  for the differences and when they matter. Ignored for other structural
  models and for `n_steps = 1`.

- n_cores:

  Positive integer. Number of processes to spread the random starts
  over. Default `1` (sequential).

- ...:

  Additional arguments passed to the measurement or structural model
  constructors (e.g., `max_val` for multinoulli models).

## Value

An object of class `mixture_model`, a list with:

- `n_components` Number of latent classes.

- `weights` Numeric vector of estimated class proportions.

- `mm` Fitted measurement model state object.

- `sm` Fitted structural model state object, or `NULL`.

- `metrics` Named list: `ll` (log-likelihood), `aic`, `bic`, `sabic`,
  `n_params`, and `entropy` (relative entropy, 0-1 scale).

- `log_resp` Matrix of log posterior class probabilities (n x K). Use
  `exp(fit$log_resp)` to obtain posterior probabilities.

- `converged` Logical. Whether the EM algorithm converged.

- `n_iter` Integer. Number of EM iterations run.

- `step1_metrics` Named list of Step-1 fit indices (only when
  `n_steps = 3`).

## Examples

``` r
# Binary LCA with 3 classes and 5 random restarts
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 3, measurement = "binary", n_init = 5)
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
print(fit)
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 3
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 252 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -337.02
#>   Parameters     : 17
#>   AIC            : 708.04
#>   BIC            : 752.33
#>   SABIC          : 698.64
#>   Rel. Entropy   : 0.4577
#>   Best solution  : found by 5 of 5 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 47.50%
#>   Class 2: 29.86%
#>   Class 3: 22.64%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
summary(fit)
#> Notice: No structural model found. Use measurement_summary() for item parameters.
measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2 | Class 3
#> -------------------------------------------------- 
#> Item_1               |   0.345 |   0.724 |   0.441
#> Item_2               |   0.737 |   0.462 |   0.229
#> Item_3               |   0.266 |   0.829 |   0.027
#> Item_4               |   0.635 |   0.366 |   0.083
#> Item_5               |   0.473 |   0.367 |   0.599
#> 
#> The Overall column, holding the observed marginal for each item, is omitted above: this fit either does not store its raw indicators - in which case refitting with the current version enables it - or holds item parameters that cannot be matched to them by name, as a multiple-group measurement model does.
#> =========================================================

# Continuous LPA (2 classes)
X_cont <- matrix(rnorm(300), nrow = 100)
fit_lpa <- fit_mixture(X_cont, n_components = 2, measurement = "continuous")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.

# 3-step LCA with a covariate and ML correction
Z <- matrix(rnorm(100), nrow = 100)
fit_cov <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
                       structural = "covariate",
                       n_steps = 3, correction = "ML", n_init = 5)
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
summary(fit_cov)
#> =========================================================
#>              STRUCTURAL MODEL SUMMARY                    
#> =========================================================
#> 
#> CATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)
#> Reference Class: 1
#> Standard errors: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                               OR         [95% CI]         P-Value
#> 
#> Class 2 ON
#>   Intercept                0.725  [    0.014,    37.579]     0.873
#>   V1                       0.780  [    0.288,     2.114]     0.625
#> =========================================================
```
