# Bootstrap Likelihood Ratio Test (BLRT) for Class Enumeration

Tests whether a model with more classes fits significantly better than
one with fewer, using a parametric bootstrap to approximate the null
distribution of the likelihood-ratio statistic. This avoids the known
violation of standard chi-squared regularity conditions in mixture
models, where the null places a parameter on the boundary of the
parameter space.

Both `blrt()` and `calc_blrt()` fit the smaller- and larger-class models
on the observed data, compute the observed likelihood-ratio statistic,
then generate `n_reps` synthetic datasets under the smaller model to
build the reference distribution. `blrt()` is the preferred name;
`calc_blrt()` is retained for backward compatibility.

## Usage

``` r
blrt(
  indicators,
  k_small,
  k_large,
  measurement = "binary",
  n_reps = 100,
  n_init_base = 20,
  n_init_boot = 10,
  verbose = TRUE,
  ...,
  X = NULL
)

calc_blrt(
  X,
  k_small,
  k_large,
  measurement = "binary",
  n_reps = 100,
  n_init_base = 20,
  n_init_boot = 10,
  verbose = TRUE,
  ...
)
```

## Arguments

- indicators:

  Matrix or data frame of measurement items. (`X` is accepted as a
  deprecated alias.)

- k_small:

  Number of classes in the smaller (null) model.

- k_large:

  Number of classes in the larger (alternative) model. Must be strictly
  greater than `k_small`.

- measurement:

  Measurement specification, as in
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  (a single type string or a named list for mixed-type indicators).

- n_reps:

  Number of bootstrap replications. Default `100`.

- n_init_base:

  Random restarts when fitting the observed-data models. Default `20`.

- n_init_boot:

  Random restarts per bootstrap replicate. Default `10`.

- verbose:

  Logical; print progress while bootstrapping. Default `TRUE`.

- ...:

  Additional arguments passed to the fitting engine.

- X:

  Deprecated alias for `indicators`.

## Value

An object of class `blrt_test`: a list with `p_value`, `obs_diff` (the
observed \\2\\\Delta\ell\\ statistic), `null_dist` (the bootstrap null
distribution), and the compared class counts. It has `print` and `plot`
methods.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
res <- blrt(X, k_small = 2, k_large = 3, measurement = "binary", n_reps = 50)
res                 # clean printed summary
plot(res)           # null distribution with the observed statistic marked
res$p_value
} # }
```
