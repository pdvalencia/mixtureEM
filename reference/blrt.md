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

A growth model is specified by a design matrix in time and a covariance
structure rather than by a measurement string, which is more than a
`measurement =` argument can carry. Pass the fitted model itself with
`from_fit =` instead: the data, the time design, the random effects and
the covariance constraints are all read off it, so the null and
alternative models differ from the fit in the number of classes and in
nothing else. That is the condition the test requires: likelihood-ratio
tests "compare models that differ only in the number of classes ... but
are not appropriate for comparing models that allow for different types
of between-class differences" (Ram & Grimm, 2009, p. 571). Passing the
fit guarantees it by construction.

## Usage

``` r
blrt(
  indicators,
  k_small,
  k_large,
  measurement,
  n_reps = 100,
  n_init_base = 20,
  n_init_boot = 10,
  verbose = TRUE,
  ...,
  from_fit = NULL,
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
  deprecated alias.) Not needed when `from_fit` is given.

- k_small:

  Number of classes in the smaller (null) model.

- k_large:

  Number of classes in the larger (alternative) model. Must be strictly
  greater than `k_small`.

- measurement:

  Measurement specification, as in
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  (a single type string or a named list for mixed-type indicators).
  Required, except when `from_fit` is supplied, which reads the
  specification off the fitted model.

- n_reps:

  Number of bootstrap replications. Default `100`, following Dziak et
  al. (2014) and the general advice in Davison and Hinkley
  (1997, p. 143) that the number of replicates be at least 99. It is a
  floor rather than a target: the attainable p-values are \\1/(B+1),
  2/(B+1), \ldots\\, so a decision that turns on the third decimal needs
  `n_reps = 999`, which recovers about 0.95 of the power of the full
  test where 99 draws recover about 0.83 (and only about 0.60 at
  \\\alpha = .01\\). See
  [`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md).

- n_init_base:

  Random restarts when fitting the observed-data models. Default `20`,
  as elsewhere in the package.

- n_init_boot:

  Random restarts per bootstrap replicate. Default `10`. This is a
  compute compromise rather than a recommended value: the two models are
  refitted `2 * n_reps` times, so the replicate search is where the cost
  of the test lives. Dziak et al. (2014) used 50 and note that too few
  restarts under the alternative can make the likelihood ratio come out
  negative. `blrt()` counts those draws and warns when there are any; if
  it does, raise this to `50`.

- verbose:

  Logical; print progress while bootstrapping. Default `TRUE`.

- ...:

  Additional arguments passed to the fitting engine.

- from_fit:

  A model fitted by
  [`fit_gmm`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  or
  [`fit_lcga`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
  whose data and specification are used for both models and every
  replicate. `k_small` and `k_large` still say which class counts to
  compare; everything else comes from the fit.

- X:

  Deprecated alias for `indicators`.

## Value

An object of class `blrt_test`: a list with `p_value`, `obs_diff` (the
observed \\2\\\Delta\ell\\ statistic), `null_dist` (the bootstrap null
distribution), `n_negative` (how many draws produced a negative
statistic), and the compared class counts. It has `print` and `plot`
methods.

## References

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their
Application* (ch. 4). Cambridge University Press.

Dziak, J. J., Lanza, S. T., & Tan, X. (2014). Effect size, statistical
power and sample size requirements for the bootstrap likelihood ratio
test in latent class analysis. *Structural Equation Modeling*, *21*(4),
534-552.
[doi:10.1080/10705511.2014.919819](https://doi.org/10.1080/10705511.2014.919819)

McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test
statistic for the number of components in a normal mixture. *Applied
Statistics*, *36*(3), 318-324.
[doi:10.2307/2347790](https://doi.org/10.2307/2347790)

Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
identifying differences in longitudinal change among unobserved groups.
*International Journal of Behavioral Development*, *33*(6), 565-576.
[doi:10.1177/0165025409343765](https://doi.org/10.1177/0165025409343765)

Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, *14*(4),
535-569.
[doi:10.1080/10705510701575396](https://doi.org/10.1080/10705510701575396)

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
