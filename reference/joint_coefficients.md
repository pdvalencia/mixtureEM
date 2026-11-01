# Association coefficients of a two-latent-variable model

The multinomial logit of `cs` on `ca` from
[`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md),
in the reference-class form: every log-odds is against the last `cs`
class, and `ca` enters as indicators against the last `ca` class. Rows
are the intercept (`term = "Intercept"`, the log-odds among the last
`ca` class), the `ca` effects (`term = "CA#k"`), and each covariate on
`cs`: one row for a common slope, or, for a covariate acting by `ca`
class, its slope in the last `ca` class (`ca = NA`) and the difference
from it in each other class.

## Usage

``` r
joint_coefficients(fit)
```

## Arguments

- fit:

  A `joint_classes` fit.

## Value

A data frame with columns `cs`, `term`, `ca`, `estimate`, `se`, `z`,
`p`. Each row is also a parameter of
[`wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_test.md)
and
[`derived_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/derived_estimates.md),
under the label `"CS#j Intercept"`, `"CS#j ON CA#k"`, `"CS#j ON x"` or,
for the difference of a covariate's slope in `ca` class k from the last
class's, `"CS#j ON x | CA#k"`.
