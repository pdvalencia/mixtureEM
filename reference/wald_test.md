# Wald Test of Linear Constraints on Model Parameters

Tests a set of linear equality constraints on the estimated parameters
of a fitted model – for example, that a growth-factor mean is the same
in three classes – without refitting the model under the constraints.
The constraints name parameters by label; the labels are those of
[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md)
for a growth model and those of
[`vcov()`](https://pdvalencia.github.io/mixtureEM/reference/vcov.mixture_model.md)
for a covariate model (see Details).

## Usage

``` r
wald_test(fit, constraints, rhs = 0, se = c("robust", "hessian", "opg"))
```

## Arguments

- fit:

  A model from
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  or
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
  or a `mixture_model` with a covariate structural model.

- constraints:

  A list or matrix of constraints. See Details.

- rhs:

  Numeric, the right-hand side \\r\\: one value per constraint, or a
  single value recycled. Default `0`.

- se:

  For growth models, which covariance estimator to use: `"robust"`
  (default), `"hessian"` or `"opg"`, as in
  [`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md).
  Ignored for covariate models, whose covariance is the one the model
  was fitted with (see
  [`vcov()`](https://pdvalencia.github.io/mixtureEM/reference/vcov.mixture_model.md)).

## Value

A one-row data frame of class `mixture_wald` with columns `Wald_Chi2`,
`df`, `p_value` and `Method`, the name of the covariance estimator
behind the statistic.

## Details

Writing the constraints as \\R\theta = r\\, the statistic is \$\$W =
(R\hat\theta - r)^\top (R \hat V R^\top)^{-1} (R\hat\theta - r),\$\$
referred to a chi-squared distribution with as many degrees of freedom
as there are constraints (Wald, 1943). \\\hat V\\ is the estimated
covariance matrix of the parameters.

**Labels.** For a
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
fit, every row of `growth_estimates(fit)` is one parameter, labelled
`"<parameter>[<item>, <class>]"`: `"growth_mean[intercept, 1]"`,
`"growth_variance[linear, 2]"`, `"residual_variance[t1, 1]"`, and so on.
A parameter constrained equal across classes appears once per class with
the same value, and equating two of its copies is not a constraint at
all. For a mixture model with a covariate structural model the labels
are the names of `vcov(fit)`, `"Class k:predictor"`, which are the
multinomial-logit coefficients relative to the class reported in its
`ref_class` attribute.

**Constraints.** Either a list, each element one constraint:

- a pair of labels, `c("a", "b")`, meaning `a - b = 0`;

- a named numeric vector of weights, `c(a = 1, b = -2)`, meaning
  `a - 2 b = rhs`;

or a numeric matrix with one row per constraint whose column names are
labels. Three classes are equal on a parameter when two pairwise
equalities hold, so `list(c("a1", "a2"), c("a2", "a3"))` has 2 degrees
of freedom; a third equality, `c("a1", "a3")`, would add nothing and is
refused as redundant.

**Which test to use.** `wald_test()` answers "are these parameters equal
(or equal to these values)?" for any parameters of the model.
[`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
answers the narrower and common question "does this covariate affect
class membership at all?", and needs only the covariate's name; on a
covariate model it is the same statistic as `wald_test()` with every one
of that covariate's coefficients set to zero.
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md)
answers the same question from bootstrapped standard errors, for samples
too small to trust the analytical ones.

Like any Wald test, this one assumes the estimates are near-normal,
which fails near the edge of the parameter space: a variance close to
zero, or a class with very few cases. A likelihood-ratio test of the
constrained model against the free one settles those cases.

## References

Wald, A. (1943). Tests of statistical hypotheses concerning several
parameters when the number of observations is large. *Transactions of
the American Mathematical Society*, *54*(3), 426-482.
[doi:10.1090/S0002-9947-1943-0012401-3](https://doi.org/10.1090/S0002-9947-1943-0012401-3)

## See also

[`growth_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/growth_estimates.md)
for the labels and estimates of a growth model,
[`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md),
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md).

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
# Do the two classes share their mean slope?
wald_test(fit, list(c("growth_mean[linear, 1]", "growth_mean[linear, 2]")))
#> =========================================================
#>                 WALD TEST (CONSTRAINTS)                  
#> =========================================================
#>   Method    : Sandwich (robust)
#> ---------------------------------------------------------
#>   Wald χ²(1) = 109.370,  p < .001
#> =========================================================
# }
```
