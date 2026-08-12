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

For the growth models the comparison should include the one-class
solution, which is the ordinary latent growth curve model: it is the
benchmark the class solutions have to beat, and reporting it is asked
for by name in the GRoLTS reporting checklist (van de Schoot et al.,
2017, item 11). It is therefore in the default `k_range` for `"gmm"` and
`"lcga"` and not for the other two, where a one-status LTA is not a
model anyone reports.

The table informs the decision; it does not make it. As Ram and Grimm
(2009, p. 571) put it, "there is not a deterministic set of rules to
follow when selecting the best model. Rather, model selection is an art
— informed by theory, past findings, past experience, and a variety of
statistical fit indices."

**Reading the `Entropy` column.** Relative entropy describes how cleanly
a solution separates its classes. The usual anchors are 0.40, 0.60 and
0.80 for low, medium and high separation (Clark & Muthen, 2009, as
reported by Lee et al., 2023, p. 653), and Ram and Grimm (p. 571)
suggest preferring the higher-entropy model when choosing among models
with similar BIC. Those anchors are on the same normalisation this
package uses. Entropy is not evidence for how many classes there are,
though, and there is no threshold it has to clear — "there are no set
cut-off criteria for deciding whether the entropy is reasonably high"
(Jung & Wickrama, 2008, p. 312) — so the package applies none.

**Reading the `Unreplicated` column.** `TRUE` means that K's reported
maximum was found by exactly one random start; refit those with
`n_init = 100` before reporting. The per-model warning is suppressed
inside this loop, since it would otherwise fire once per K.

## Usage

``` r
compare_longitudinal(
  indicators,
  k_range = NULL,
  model = c("lta", "rmlca", "gmm", "lcga"),
  times = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- indicators:

  The repeated indicators, in any format accepted by
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md);
  for `model = "gmm"` or `"lcga"`, the single repeated outcome, in any
  format accepted by
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md).

- k_range:

  Integer vector of class or status counts to fit. Defaults to `1:4` for
  the growth models and `2:4` for the others.

- model:

  `"lta"` (default), `"rmlca"`, `"gmm"` or `"lcga"`.

- times:

  Number of occasions; required for wide input.

- verbose:

  Print progress.

- ...:

  Further arguments passed to the fitting function, such as
  `measurement`, `time_invariance` or `tau_homogeneous` for the
  categorical models, `degree`, `time_scores`, `random_effects`, `psi`
  or `residual` for the growth models, and `n_init` for any of them.
  Passing the same specification to every K is the point: the criteria
  in the table are only comparable across models that differ in nothing
  else. Ram and Grimm (2009, p. 571) state the same rule for the tests:
  they "compare models that differ only in the number of classes ... but
  are not appropriate for comparing models that allow for different
  types of between-class differences".

## Value

A list with `fit_table` (columns `Classes`, `LL`, `Params`, `AIC`,
`BIC`, `SABIC`, `Entropy` and `Unreplicated`), the fitted `models`
(named `"K2"`, `"K3"`, ...) and `best_k`, the class count with the
lowest BIC.

## References

van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., &
Vermunt, J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on
latent trajectory studies. *Structural Equation Modeling*, *24*(3),
451-467.

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
