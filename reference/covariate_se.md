# Standard Errors for Covariate Effects in Three-Step Models

Relating latent class membership to covariates in two or three steps
makes the third step's "data" — the posterior class assignments, and
under the ML adjustment the classification table built from them —
estimates rather than observations. Ignoring that gives confidence
intervals which are too narrow. The `se` argument of
[`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
selects among three estimators of the covariate covariance matrix.

## Details

- `"corrected"` (default):

  The first-order corrected estimator of Bakk et al. (2014, eq. 17),
  \\D_3^\* = D_3 + J D_1 J'\\, where \\D_3\\ is the step-3 sandwich
  below, \\D_1\\ the sampling variance of the step-1 measurement model,
  and \\J = (-H_3)^{-1} \partial^2 L_3 / \partial\theta_3
  \partial\theta_1'\\ the implicit derivative of the step-3 estimates
  with respect to the step-1 parameters. This is the estimator their
  simulation recommends; the second-order form of their eq. 18 adds a
  term that vanishes asymptotically. \\D_1\\ needs the step-1 model in
  an unconstrained parameterisation, which exists for binary,
  categorical, count and continuous indicators, alone, mixed, in time
  blocks, by group, or with item DIF. For a growth mixture model or an
  LCGA as step 1, the `"robust"` estimator is reported instead, with a
  warning, and [`summary()`](https://rdrr.io/r/base/summary.html) names
  it.

- `"robust"`:

  The step-3 sandwich only, \\D_3 = (-H_3)^{-1} M (-H_3)^{-1}\\, with
  \\H_3\\ the marginal step-3 Hessian and \\M\\ the outer product of the
  case-level scores (PSU-level within strata when a survey design is
  attached). Proportional assignment gives each case \\K\\ weighted
  records, so this sandwich is needed even before any step-1 uncertainty
  is considered.

- `"hessian"`:

  \\(-H_3)^{-1}\\ alone. Provided for completeness; it ignores the
  record duplication and is not recommended.

None of these is the Hessian the package reported previously, which came
from the M-step and measured the curvature of the *Q function* rather
than of the step-3 log-likelihood, and so was smaller still.

`"corrected"` is the statistically right answer and is the default for
that reason. `"robust"` is the step-3 sandwich alone, the estimator of
an analysis that treats the assigned classes as given; use
`se = "robust"` when that is the analysis being reported. Under the
default the intervals are wider, and the difference is a difference in
estimator — the corrected form carries step-1 uncertainty that the
sandwich alone omits.

How much this matters depends almost entirely on how well separated the
classes are. A 250-replication coverage study on the design of Bakk et
al. (`data-raw/covariate_se_simulation.R` in the package sources) gives,
for a nominal 95 percent interval, ranged over six covariate
coefficients:

|                     |                |               |            |               |
|---------------------|----------------|---------------|------------|---------------|
| **design**          | **entropy R2** | **Q Hessian** | **robust** | **corrected** |
| n = 500, rho = .80  | .64            | .54 - .88     | .83 - .98  | .92 - .98     |
| n = 500, rho = .90  | .88            | .87 - .94     | .94 - .98  | .94 - .98     |
| n = 2000, rho = .90 | .88            | .87 - .95     | .94 - .97  | .94 - .97     |

where the first column is the estimator this package reported before the
three above existed. With poorly separated classes the correction is the
difference between an interval that covers and one that does not, and it
bites hardest on the largest effects. Once entropy R-squared reaches
about .90 the corrected and robust estimators agree to the third
decimal, which matches Bakk et al.'s finding that above n = 2000 with
entropy over .90 no correction is needed.

A caveat the same study makes visible: at low separation the three-step
*point* estimates are themselves biased (Bakk et al., discussion),
because step 1 underestimates the classification error. No variance
estimator repairs that. A covariate coefficient from a model with
entropy below about .60 should be treated with suspicion however wide
its interval.

[`summary()`](https://rdrr.io/r/base/summary.html),
[`confint.mixture_model`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md)
and
[`analytical_wald_test`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
all print the name of the estimator that produced the numbers they show.

## The two-step estimator

A two-step fit (`n_steps = 2`, or `steps = 2` in
[`add_covariates`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md))
has no assigned classes and no classification table: the
class-prediction coefficients maximise the joint log-likelihood with the
measurement parameters held at their step-1 estimate. The same three
`se` values select among the pseudo-maximum-likelihood estimators of
Gong and Samaniego (1981) as written for this case by Bakk and Kuha
(2018, eq. 5):

- `"corrected"` (default):

  \\V = V_2 + V_1\\, where \\V_2 = I\_{22}^{-1}\\ is the inverse
  observed information of the joint log-likelihood in the structural
  parameters, and \\V_1 = I\_{22}^{-1} I\_{12}' \Sigma\_{11} I\_{12}
  I\_{22}^{-1}\\ carries the sampling variance \\\Sigma\_{11}\\ of the
  step-1 measurement parameters through the cross-information
  \\I\_{12}\\. Both blocks of the information matrix are differenced
  numerically from the joint log-likelihood; \\\Sigma\_{11}\\ is the
  step-1 observed information inverted, or the outer product of the
  step-1 scores above 100 measurement parameters.

- `"robust"`:

  The sandwich \\V_2 M V_2\\ with \\M\\ the outer product of the
  case-level scores of the joint log-likelihood (PSU-level within strata
  under a survey design), without the step-1 term.

- `"hessian"`:

  \\V_2\\ alone: what an analysis that treats the step-1 estimate as
  known would report.

The step-1 term is not small when the classes are poorly separated. Bakk
and Kuha (2018, Table 2) put the step-2-only standard error at 40 to 60
percent of the whole at entropy R-squared near 0.36, 73 to 88 percent
near 0.65 and 98 to 99 percent near 0.90, with the 95 percent intervals
of the step-2-only estimator covering 77 to 86 percent at the lowest
separation. Under a survey design `"corrected"` uses the model-based
\\V_2\\, not the linearised one; `"robust"` is the design-based choice
there.

## Scope

The corrected and robust estimators cover a covariate (class-prediction)
structural model estimated with `n_steps = 3` and `correction = "none"`
or `"ML"`, and with `n_steps = 2`. Four cases fall back to the
uncorrected Hessian, and say so in the printed output:
`correction = "BCH"` (whose weights need their own variance treatment,
and which is not recommended for covariates in any case); `n_steps = 1`,
where measurement and structural parameters are estimated jointly and no
carry-over correction applies; a two-step fit whose structural model is
a distal outcome rather than class predictors; and a covariate combined
with a distal outcome in one nested structural model. The step-1 term
additionally requires a measurement model whose parameters this package
can put on an unconstrained scale — binary, polytomous, Gaussian, count,
mixed, and repeated-measures models qualify; growth models do not, and
there the robust sandwich (three-step) or \\V_2\\ (two-step) is reported
instead.

## References

Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
assignments to external variables: Standard errors for correct
inference. *Political Analysis*, *22*(4), 520-540.
[doi:10.1093/pan/mpu003](https://doi.org/10.1093/pan/mpu003)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450-469.

Gong, G., & Samaniego, F. J. (1981). Pseudo maximum likelihood
estimation: Theory and applications. *The Annals of Statistics*, *9*(4),
861-869.

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, *83*(4),
871-892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

## See also

[`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
[`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md),
[`confint.mixture_model`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md),
[`analytical_wald_test`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
