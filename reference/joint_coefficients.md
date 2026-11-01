# Association coefficients of a two-latent-variable model

The multinomial logit of `y` on `x` from
[`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md),
in the reference-class form: every log-odds is against the last `y`
class, and `x` enters as indicators against the last `x` class. Rows are
the intercept (`term = "Intercept"`, the log-odds among the last `x`
class), the `x` effects (`term = "X#k"`), and each covariate on `y`: one
row for a common slope, or, for a covariate acting by `x` class, its
slope in the last `x` class (`x = NA`) and the difference from it in
each other class.

## Usage

``` r
joint_coefficients(fit)
```

## Arguments

- fit:

  A `joint_classes` fit.

## Value

A data frame with columns `y`, `term`, `x`, `estimate`, `se`, `z`, `p`.
Each row is also a parameter of
[`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md)
and
[`derived_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/derived_estimates.md),
under the label `"Y#j Intercept"`, `"Y#j ON X#k"`, `"Y#j ON x"` or, for
the difference of a covariate's slope in `x` class k from the last
class's, `"Y#j ON x | X#k"`.
