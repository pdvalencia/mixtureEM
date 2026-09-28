# Direct covariate effects on the indicators

Reports the slopes a model fitted with `predictors_items` estimates: the
effect of a covariate on the log-odds of endorsing an item, for people
in the same latent class. A non-zero slope is differential item
functioning (DIF): the item does not measure the classes the same way
for everyone. A *uniform* slope is one number shared by every class; a
*class-specific* slope has one value per class.

The item probabilities the fit reports
([`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)) are those of a
case whose item covariates are all zero, so centre or code a covariate
so that zero is a meaningful reference.

## Usage

``` r
dif_effects(fit, se = TRUE)
```

## Arguments

- fit:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  with `predictors_items`.

- se:

  Logical; compute standard errors (costs one numerical Hessian of the
  full parameter vector).

## Value

A data frame with one row per slope: `item`, `covariate`, `class`
(`"all"` for a uniform slope), `estimate`, `se`, `z` and `p`.

## Details

Which pairs to include is the analyst's choice; the package does not
search for them. The usual workflow (Masyn, 2017) fits the model without
direct effects, inspects the covariate-item residuals, adds the pairs
they point to, and compares the nested fits with
[`lr_test`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).
Two cautions from the literature: a DIF screen produces false positives
when local independence is also violated, increasingly so at large
samples (Depaoli, Jia & Visser, 2025), and a class-specific slope needs
considerably larger samples and classes to detect than a uniform one,
which is why uniform is the default.

Standard errors invert the numerical observed information of the full
one-step log-likelihood: item thresholds, direct effects and the
class-membership regression together.

## References

Depaoli, S., Jia, F., & Visser, I. (2025). *Structural Equation
Modeling*, 32(5), 780-800.

Masyn, K. E. (2017). Measurement invariance and differential item
functioning in latent class analysis with stepwise multiple indicator
multiple cause modeling. *Structural Equation Modeling*, 24(2), 180-197.

## See also

[`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
[`lr_test`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
