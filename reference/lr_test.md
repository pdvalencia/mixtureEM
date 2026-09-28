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

- **The same question before any transitions are modelled**: the
  configural against the invariant measurement model, with the occasions
  unlinked, which is how a three-step analysis settles invariance before
  step 3 (Nylund-Gibson, Arch & Carter, 2026). A three-step fit's
  `$step1` is exactly that model, so fit
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `n_steps = 3` and `measurement_invariance = "full"` and `"none"`,
  and test `lr_test(fit_full$step1, fit_none$step1)`. With `"none"` the
  step-1 model is one free latent class model per occasion; with
  `"full"` the item probabilities are shared, so the degrees of freedom
  are the item parameters freed at every occasion after the first.

- **A time-homogeneous transition matrix** (sec. 7.14): fit
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `transition_invariance = "full"` and `"none"`.

The models must be nested and fitted to the same data. That is not
checked beyond the parameter counts and sample size, so it remains the
analyst's responsibility.

A three-step fit (`fit_lta(n_steps = 3)`) reports the log-likelihood of
its third step, which models the assigned statuses rather than the
items, so it is refused beside a one- or two-step fit. Two three-step
fits are tested only when their third steps saw the same data: the same
assigned statuses and the same classification-error matrices.

Because `full` strictly nests `restricted`, its log-likelihood can never
be genuinely lower — if it comes out that way here, the `full` model's
random-restart search landed on a worse local optimum than the
`restricted` model's did, not a real result. A warning is issued in that
case; refitting `full` with a larger `n_init` is the usual fix.

## Usage

``` r
lr_test(restricted, full, scaled = c("auto", "yes", "no"))
```

## Arguments

- restricted:

  The more constrained model (fewer parameters).

- full:

  The less constrained model.

- scaled:

  Whether to apply the Satorra-Bentler/Asparouhov scaling correction.
  `"auto"` (the default) applies it exactly where it is required for
  validity — under sampling weights or a complex survey design — and is
  the historical behaviour. `"yes"` applies it to an unweighted pair as
  well, giving the robust (MLR-scaled) difference test; `"no"`
  suppresses it even under weights. A negative scaled statistic means
  the correction has failed for that pair, and `statistic_raw` should be
  reported instead.

## Value

A list of class `"lr_test"`.

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent Class and Latent
Transition Analysis: With Applications in the Social, Behavioral, and
Health Sciences*. Wiley.

Nylund-Gibson, K., Arch, D. A. N., & Carter, D. (2026). Latent
transition analysis with auxiliary variables: A demonstration of the ML
3-step and BCH in Mplus. *The Quantitative Methods for Psychology*,
*22*(1).
[doi:10.20982/tqmp.22.1.p001](https://doi.org/10.20982/tqmp.22.1.p001)
