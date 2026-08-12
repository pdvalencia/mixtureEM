# Compare Mixture Models Across a Range of Class Numbers

Fits a sequence of measurement-only mixture models, one for each value
of `k` in `k_range`, and returns a table of fit indices to guide class
enumeration. The best model according to BIC is identified
automatically.

## Usage

``` r
compare_mixtures(
  X,
  k_range = 1:5,
  measurement = "binary",
  n_init = 10,
  n_steps = 1,
  ...
)
```

## Arguments

- X:

  A numeric matrix or data frame of indicator variables.

- k_range:

  Integer vector of class numbers to fit. All values must be \>= 1.
  Default is `1:5`.

- measurement:

  Character string specifying the measurement model type. See
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  for accepted values. Default is `"binary"`.

- n_init:

  Positive integer. Number of random restarts per model. Default is
  `10`.

- n_steps:

  Integer. Estimation method: `1`, `2`, or `3`. Default is `1`.

- ...:

  Additional arguments passed to
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

## Value

A named list with three elements:

- `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
  `Params`, `AIC`, `BIC`, `SABIC`, `Entropy` and `Unreplicated`.

- `models` Named list of fitted `mixture_model` objects, one per K
  (names are `"K1"`, `"K2"`, etc.).

- `best_k` Integer. The value of K with the lowest BIC.

## Details

**Reading the `Entropy` column.** Relative entropy describes how cleanly
a solution separates the classes, on a 0-to-1 scale. The usual anchors
are 0.40, 0.60 and 0.80 for low, medium and high separation (Clark &
Muthen, 2009, as reported by Lee et al., 2023, p. 653); Ram and Grimm
(2009, p. 571) put the same point as "high values of entropy (\>.80)
indicate that individuals are classified with confidence", and suggest
preferring the higher-entropy model when choosing among models with
similar BIC. Those anchors are on the same normalisation this package
uses.

Entropy is not evidence for how many classes there are, and there is no
threshold it has to clear: "there are no set cut-off criteria for
deciding whether the entropy is reasonably high" (Jung & Wickrama, 2008,
p. 312). The numbers above are for reading a table, and mixtureEM
applies no entropy threshold anywhere.

**Reading the `Unreplicated` column.** `TRUE` means the reported maximum
for that K was found by exactly one random start. Refit those values of
K with `n_init = 100` before reporting them; see
[`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md).

## References

Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the
number of classes in latent class analysis and growth mixture modeling:
A Monte Carlo simulation study. *Structural Equation Modeling*, *14*(4),
535-569.
[doi:10.1080/10705510701575396](https://doi.org/10.1080/10705510701575396)

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

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
                           n_init = 5)
#> Running Model Selection across K = 1 to 4...
#> 
#> Fitting 1-class model...
#> Fitting 2-class model...
#> Fitting 3-class model...
#> Fitting 4-class model...
#> 
#> === Model Selection Summary ===
#>   Classes       LL Params     AIC     BIC   SABIC Entropy Unreplicated
#> 1       1 -342.102      5 694.203 707.229 691.438   1.000        FALSE
#> 2       2 -340.085     11 702.169 730.826 696.085   0.303        FALSE
#> 3       3 -337.020     17 708.040 752.328 698.637   0.458        FALSE
#> 4       4 -334.543     23 715.086 775.005 702.365   0.590        FALSE
#> 
#> -> Best model according to BIC: 1 classes
result$fit_table
#>   Classes        LL Params      AIC      BIC    SABIC   Entropy Unreplicated
#> 1       1 -342.1016      5 694.2032 707.2290 691.4378 1.0000000        FALSE
#> 2       2 -340.0846     11 702.1692 730.8261 696.0854 0.3033780        FALSE
#> 3       3 -337.0199     17 708.0398 752.3277 698.6374 0.4575710        FALSE
#> 4       4 -334.5428     23 715.0857 775.0046 702.3648 0.5898679        FALSE
result$best_k
#> [1] 1
```
