# Derived Parameters with Standard Errors and Intervals

Estimates functions of a fitted model's parameters – a sum of two
coefficients, a difference, a probability built from several of them –
each with a standard error, a z test against zero and a confidence
interval. The parameters are named by the labels
[`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md)
uses.

## Usage

``` r
derived_estimates(fit, ..., level = 0.95, se = c("robust", "hessian", "opg"))
```

## Arguments

- fit:

  A model
  [`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md)
  accepts.

- ...:

  Named derived parameters: weight vectors or functions (see Details).

- level:

  Confidence level of the intervals.

- se:

  For growth models, which covariance estimator to use, as in
  [`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md).
  Ignored otherwise.

## Value

A data frame with one row per derived parameter: `parameter`,
`estimate`, `se`, `z`, `p`, `lower`, `upper`.

## Details

Each derived parameter is given as one named argument in `...`, either

- a named numeric vector of weights, `c(a = 1, b = 1, c = -1)`, meaning
  \\a + b - c\\; its variance is \\c^\top V c\\, exact; or

- a function of the named vector of estimates, returning one number,
  such as `function(p) plogis(p[["a"]] + p[["b"]])`. Its standard error
  is the delta method's, \\\sqrt{g^\top V g}\\ with \\g\\ the gradient
  of the function at the estimates, taken by central differences.

\\V\\ is the covariance matrix the fit reports (for a growth model, the
one `se` selects). The interval is the estimate plus or minus the normal
quantile times the standard error, on the scale of the derived
parameter. For a probability near 0 or 1 an interval built on the logit
scale and transformed back stays inside the unit interval; derive the
logit, and transform its bounds with
[`plogis()`](https://rdrr.io/r/stats/Logistic.html).

For
[`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md)
the labels are those of
[`joint_coefficients()`](https://pdvalencia.github.io/mixtureEM/reference/joint_coefficients.md):
`"CS#j Intercept"`, `"CS#j ON CA#k"`, `"CS#j ON x"` and
`"CS#j ON x | CA#k"`. A covariate `x` acting on `cs` by `ca` class has
slope `"CS#j ON x"` in the last `ca` class and that plus
`"CS#j ON x | CA#k"` in class k, so the association among cases with
`x = 1` is `c("CS#j ON CA#k" = 1, "CS#j ON x | CA#k" = 1)`.
[`joint_probabilities()`](https://pdvalencia.github.io/mixtureEM/reference/joint_probabilities.md)
builds the conditional class probabilities this way.

## References

Oehlert, G. W. (1992). A note on the delta method. *The American
Statistician*, *46*(1), 27-29.

## See also

[`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md)
to test several restrictions at once,
[`joint_probabilities()`](https://pdvalencia.github.io/mixtureEM/reference/joint_probabilities.md).

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 800
w  <- rbinom(n, 1, 0.5)
ca <- sample(1:2, n, TRUE)
cs <- ifelse(runif(n) < plogis(-1 + 2 * (ca == 1) + w * (ca == 1)), 1L, 2L)
A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(ca == 1, 0.85, 0.15)))
S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cs == 1, 0.85, 0.15)))
fit <- fit_joint_classes(A, S, predictors_cs = data.frame(w = w),
                         cs_effects = "by_ca", n_init = 5, random_state = 1)
# The association among cases with w = 1, and the probability of cs class 1
# for ca class 1 at w = 1.
derived_estimates(fit,
  assoc_w1 = c("CS#1 ON CA#1" = 1, "CS#1 ON w | CA#1" = 1),
  p_w1 = function(p) plogis(p[["CS#1 Intercept"]] + p[["CS#1 ON CA#1"]] +
                            p[["CS#1 ON w"]] + p[["CS#1 ON w | CA#1"]]))
#>   parameter  estimate        se         z            p     lower    upper
#> 1  assoc_w1 2.9233534 0.3441608  8.494149 1.993885e-17 2.2488106 3.597896
#> 2      p_w1 0.6976749 0.0393972 17.708741 3.590148e-70 0.6204578 0.774892
# }
```
