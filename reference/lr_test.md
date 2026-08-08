# Likelihood-Ratio Test for Two Nested Models

Compares two nested models by the likelihood-ratio difference test,
\\-2(\ell_0 - \ell_1)\\ on \\P_1 - P_0\\ degrees of freedom. It accepts
any pair of fits from
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
or
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
and answers questions of the form "does freeing these parameters buy a
significantly better fit?"

- **Measurement invariance across groups** (Collins & Lanza, 2010, sec.
  5.8): fit
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  with `group_effects = "prevalence"` and with `"both"`, and compare.
  This is a cross-sectional test.

- **Equal prevalences across groups** (sec. 5.11): compare
  `group_effects = "none"` against `"prevalence"`.

- **Measurement invariance across time** (sec. 7.11): fit
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `measurement_invariance = "full"` and `"none"` and compare.

- **A time-homogeneous transition matrix** (sec. 7.14): fit
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `transition_invariance = "full"` and `"none"`.

The models must be nested and fitted to the same data. That is not
checked beyond the parameter counts and sample size, so it remains the
analyst's responsibility.

Because `full` strictly nests `restricted`, its log-likelihood can never
be genuinely lower — if it comes out that way here, the `full` model's
random-restart search landed on a worse local optimum than the
`restricted` model's did, not a real result. A warning is issued in that
case; refitting `full` with a larger `n_init` is the usual fix.

## Usage

``` r
lr_test(restricted, full)
```

## Arguments

- restricted:

  The more constrained model (fewer parameters).

- full:

  The less constrained model.

## Value

A list of class `"lr_test"`.

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent Class and Latent
Transition Analysis: With Applications in the Social, Behavioral, and
Health Sciences*. Wiley.
