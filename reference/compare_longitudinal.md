# Compare Longitudinal Mixture Models Across a Range of Class Counts

Fits a series of models with increasing numbers of latent classes
(RMLCA) or latent statuses (LTA) and tabulates the usual selection
criteria, mirroring
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
for the cross-sectional case. Lower AIC, BIC and SABIC are better;
entropy summarises how cleanly cases are classified.

Selecting the number of latent statuses for an LTA should use the data
from every occasion at once rather than a separate cross-sectional
analysis per occasion: pooling the repeated measures gives the model
more information, so a solution can be supported longitudinally that no
single occasion would support on its own (Collins & Lanza, 2010, sec.
7.3.3).

Note that likelihood-based tests of \\K\\ against \\K-1\\ classes, such
as the bootstrap likelihood-ratio test, do not have their usual
reference distribution here, which is why only information criteria are
reported.

## Usage

``` r
compare_longitudinal(
  indicators,
  k_range = 2:4,
  model = c("lta", "rmlca"),
  times = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- indicators:

  The repeated indicators, in any format accepted by
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).

- k_range:

  Integer vector of class or status counts to fit.

- model:

  `"lta"` (default) or `"rmlca"`.

- times:

  Number of occasions; required for wide input.

- verbose:

  Print progress.

- ...:

  Further arguments passed to
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  or
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md),
  such as `measurement`, `time_invariance`, `tau_homogeneous` or
  `n_init`.

## Value

A list with `fit_table`, the fitted `models` (named `"K2"`, `"K3"`, ...)
and `best_k`, the class count with the lowest BIC.
