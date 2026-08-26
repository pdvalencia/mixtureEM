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
  measurement,
  n_init = 10,
  n_steps = 1,
  vlmr = c("none", "standard", "robust", "both"),
  n_cores = 1L,
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

  Character string or named list specifying the measurement model type.
  Required; see
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  for the accepted values and for why there is no default.

- n_init:

  Positive integer. Number of random restarts per model. Default is
  `10`.

- n_steps:

  Integer. Estimation method: `1`, `2`, or `3`. Default is `1`.

- vlmr:

  Character string. Whether to add the Vuong-Lo-Mendell-Rubin test of K
  against K+1 classes, and in which form: `"none"` (the default),
  `"standard"` for Vuong's own formulae on the ordinary covariance
  matrix, `"robust"` for the variant that substitutes the sandwich
  covariance, or `"both"`. See Details for why it is off by default.

- n_cores:

  Positive integer. Number of processes to spread the random starts
  over, within each K. Default `1` (sequential).

- ...:

  Additional arguments passed to
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

## Value

An object of class `mixture_comparison`: a named list with three
elements, which can be indexed exactly as a plain list.

- `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
  `Params`, `AIC`, `BIC`, `SABIC`, `Entropy` and `Unreplicated`. With
  `vlmr` set it also carries `VLMR_LR` and one p-value column per
  requested form (`VLMR_p`, `VLMR_p_robust`); each row tests its own K
  against the next one in the table, so the last row is `NA`.

- `models` Named list of fitted `mixture_model` objects, one per K
  (names are `"K1"`, `"K2"`, etc.).

- `best_k` Integer. The value of K with the lowest BIC.

- `vlmr` Present only when `vlmr` is set: one entry per row of the table
  holding the likelihood-ratio statistic and, for each requested form,
  the mean and standard deviation of the reference distribution
  alongside the p-value. Those two moments are what say which
  distribution produced a given p-value, and are not printed.

[`plot()`](https://pdvalencia.github.io/mixtureEM/reference/plot.mixture_comparison.md)
draws the criteria against K.

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
K with more starts (`n_init = 100` is the usual next step) before
reporting them; see
[`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md).

**Reading the `VLMR` columns.** They appear only when `vlmr` is set, and
are off by default for two reasons. One is cost: the test needs a
numerical Hessian for each model, which is quadratic in the number of
parameters, and a function people call casually should not pay that
unasked. The other is that the test does not deserve to be printed as a
matter of course. Vermunt (2024) concludes that "neither of the two
implementations yield uniformly distributed p-values under the correct
null hypothesis, indicating this test is not the best model selection
tool in mixture modeling".

The two implementations differ only in which covariance matrix of the
parameters enters the reference distribution: one program uses the
ordinary one, another the robust (sandwich) one (Vermunt, 2024). The
difference is not cosmetic: on the same data the two can return p = .00
and p = .15. The robust version's reference distribution is much more
sensitive to the particular sample, especially when the classes are
poorly separated. Neither version's p-values are uniform under the null,
so treat a VLMR result as one input among several and prefer
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
where it is affordable.

The reason is known and is not a numerical artefact. The reference
distribution is derived from a theorem (Vuong, 1989, Theorem 3.3) that
requires the parameters of the larger model to be identified at the
point where it reduces to the smaller one. A mixture never satisfies
this: the larger model reproduces the smaller only by emptying a class
or by duplicating one, and in each case some parameters vanish from the
likelihood and the information matrix is singular (Jeffries, 2003). The
test is therefore best read as a descriptive comparison rather than a
calibrated p-value.

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

Vermunt, J. K. (2024). The Vuong-Lo-Mendell-Rubin test for latent class
and latent profile analysis. *Methodology*, *20*(1), e12467.
[doi:10.5964/meth.12467](https://doi.org/10.5964/meth.12467)

Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica*, *57*(2), 307-333.

Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
components in a normal mixture. *Biometrika*, *88*(3), 767-778.

Jeffries, N. O. (2003). A note on "Testing the number of components in a
normal mixture". *Biometrika*, *90*(4), 991-994.

Imhof, J. P. (1961). Computing the distribution of quadratic forms in
normal variables. *Biometrika*, *48*(3/4), 419-426.

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
