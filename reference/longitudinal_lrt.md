# Likelihood-Ratio Test for Nested Longitudinal Models

Compares two nested models by the likelihood-ratio difference test,
\\-2(\ell_0 - \ell_1)\\ on \\P_1 - P_0\\ degrees of freedom. This is the
test Collins and Lanza use throughout chapters 7 and 8, and it answers
the questions those chapters pose:

- **Measurement invariance across time** (sec. 7.11): fit with
  `time_invariance = "full"` and `"none"` and compare.

- **Invariance of the transition matrix across time** (sec. 7.14): fit
  with `tau_homogeneous = TRUE` and `FALSE` and compare.

- **Group differences in prevalences or transitions** (sec. 8.7, 8.8):
  compare a model with the relevant parameters constrained equal across
  groups against one that frees them.

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
longitudinal_lrt(restricted, full)
```

## Arguments

- restricted:

  The more constrained model (fewer parameters).

- full:

  The less constrained model.

## Value

A list of class `"longitudinal_lrt"`.
