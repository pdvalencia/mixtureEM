# Bivariate Residuals

A local measure of fit: for each pair of categorical indicators, the
Pearson chi-square of the observed two-way table against the table the
model implies, divided by its degrees of freedom \\(R_a - 1)(R_b - 1)\\.
If the model were true, a bivariate residual should not be substantially
larger than 1.

Where
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
says whether the model fits, this says *where* it fails: a large value
flags the specific pair of items whose association the latent classes do
not reproduce, which is the conditional-independence assumption showing
its seams. This is the classic local-dependence diagnostic (Oberski et
al., 2013), and what the usual referee question about local dependence
asks for.

Unlike the absolute-fit statistics, bivariate residuals do not require
the full response-pattern table and so remain usable with many
indicators: the model-implied two-way margin follows in closed form from
conditional independence, \\P(y_a = r, y_b = s) = \sum_k \gamma_k\\
p_a(r\|k)\\ p_b(s\|k)\\.

## Usage

``` r
bivariate_residuals(
  object,
  n_reps = 0,
  n_init_boot = 10,
  n_cores = 1L,
  verbose = FALSE
)
```

## Arguments

- object:

  A model fitted by
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  or
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  with categorical indicators, or with a plain continuous measurement
  model.

- n_reps:

  Number of parametric bootstrap replicates used to calibrate the
  residuals. The default, `0`, computes the residuals alone and is the
  cheaper, uncalibrated diagnostic. `100` is the recommended working
  value and gives a Monte Carlo standard error of about 0.022 at \\p =
  0.05\\, which is adequate for flagging a pair; `500`, the number
  Oberski et al. used, is the publication-grade setting. Ignored for a
  continuous measurement model.

- n_init_boot:

  Random starts per bootstrap replicate. Replicates are refit without
  the final refinement step, since a replicate needs a residual rather
  than polished estimates.

- verbose:

  Report bootstrap progress.

## Value

For categorical indicators, an object of class `bivariate_residuals`: a
lower-triangular indicator-by-indicator matrix, `NA` on and above the
diagonal. When `n_reps > 0` a matrix of bootstrap p-values, laid out the
same way, is attached as the `"p"` attribute and printed beside each
residual. The sum of all pairwise residuals is attached as the `"total"`
attribute and printed as "Total BVR", a single headline number for how
much local dependence the model as a whole is carrying. For a continuous
measurement model, an object of class `bivariate_residuals_gaussian`
holding the modification index and expected parameter change per pair
per class, the model-implied residual covariance and correlation, and a
count of pairs where the information matrix needed the outer-product
fallback. `NULL` (with a message) when neither applies.

## How much to trust it

The chi-square reference for this statistic does not work, and the
evidence is blunt. Over the eight null conditions of Oberski et al.
(2013, Table 1) – loadings of .5 and .8 crossed with n of 200, 500, 1000
and 5000, 200 samples each – a bivariate residual referred to chi-square
rejected at a nominal 5 percent level in **zero of 200 samples in every
one of the eight**. Its empirical mean ran between 0.25 and 0.36 against
the 1 that a \\\chi^2_1\\ has, and its variance between 0.1 and 0.2
against 2. Three consequences follow, and all three matter more than the
usual hedging suggests:

- **A low bivariate residual is not evidence of good fit.** A statistic
  that never rejects when the model is true also has nothing to say when
  it is. This is the authors' own closing point.

- **The ranking is weaker than it looks.** In their Figure 1 the naive
  bivariate residual has uniformly the lowest power of the three methods
  compared. It is adequate only against a residual correlation of about
  -0.4; it needs n of 5000 or more for correlations of \\\pm 0.2\\ and
  -0.2, and it almost never detects \\\pm 0.05\\. All three methods lose
  power as the loadings grow, so well-separated classes hide local
  dependence rather than expose it.

- **With missing data a large value is ambiguous.** Under MAR the
  observed side of any residual statistic carries selection bias
  (Asparouhov & Muthen, 2015), so a large residual may be reporting the
  missingness mechanism rather than local dependence.

`n_reps` replaces the broken reference distribution with a parametric
bootstrap, which in the same simulation held between 0.020 and 0.085
against a nominal 0.05. Use it before drawing any conclusion from the
size of a residual. Each bootstrap replicate is blanked out in the same
cells as the real data, so a p-value from missing data is exact when the
data are missing completely at random and an approximation otherwise –
the replicate's missingness cannot reproduce a dependence on class or
response that the real gaps might carry.

Note that the statistic bootstrapped is the one this function computes,
\\\chi^2\\ divided by its degrees of freedom, and not a raw Pearson
\\\chi^2\\. For binary items the degrees of freedom are 1 and the two
coincide; for polytomous items they do not, so do not compare the number
printed here against a raw chi-square table. The bootstrap is applied to
whatever statistic is computed, so it is calibrated either way.

With missing data each pair is computed on the cases observing both
items, and the expected counts are evaluated at those cases' own
posterior class membership rather than at the class proportions for the
whole sample. This matters whenever the two subsets differ – an item
that is missing more often for one class than another otherwise makes an
unrelated pair of items look locally dependent, when what actually
happened is that missingness changed who is left in the comparison. This
is a pairwise-complete statistic rather than a full-information one, so
still read it as descriptive when missingness is heavy. For categorical
indicators this is the same statistic, with the same divisor, that
another program reports as a bivariate residual, and the values agree
closely on the same fit.

For a plain continuous (Gaussian) measurement model with no missing
data, a different statistic is returned instead: for each item pair and
class, the modification index (Sorbom, 1989) for freeing the
within-class residual covariance of that pair, with its variance
adjusted for the model's other parameters, following Oberski, van
Kollenburg and Vermunt (2013). It is computed separately in each class,
because a residual dependence can run in opposite directions in
different classes and a pooled statistic would average it away. Another
program's "bivariate residual" for continuous indicators is deliberately
*not* adjusted for the model's other parameters, so the two do not agree
numerically; the modification index is preferred here on the evidence of
Oberski et al.'s simulation, in which a bivariate residual referred to
chi-square gave below-nominal size and inadequate power, while the
modification index reproduced its nominal distribution and was the more
powerful of the two adequate methods. Detection is reliable only when
the offending effects are few and the measurement model is strong
(Janssen, van Laar, de Rooij, Kuha & Bakk, 2019); unmodelled
within-class dependence is worth taking seriously mainly because it can
manufacture spurious classes (Bauer & Curran, 2004).

## References

Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
Carlo evaluation of three methods to detect local dependence in binary
data latent class models. *Advances in Data Analysis and
Classification*, *7*(3), 267-279.
[doi:10.1007/s11634-013-0146-2](https://doi.org/10.1007/s11634-013-0146-2)

Sorbom, D. (1989). Model modification. *Psychometrika*, *54*, 371-384.

Janssen, J. H. M., van Laar, S., de Rooij, M., Kuha, J., & Bakk, Z.
(2019). The detection of local dependence in the presence of one
continuous latent variable: A comparison of different statistics.
*Structural Equation Modeling*, *26*(2), 280-290.

Bauer, D. J., & Curran, P. J. (2004). The integration of continuous and
discrete latent variable models: Potential problems and promising
opportunities. *Psychological Methods*, *9*(1), 3-29.

van Kollenburg, G. H., Mulder, J., & Vermunt, J. K. (2015). Assessing
model fit in latent class analysis when asymptotics do not hold.
*Methodology*, *11*(2), 65-79.

Asparouhov, T., & Muthen, B. (2015). Residual associations in latent
class and latent transition analysis. *Structural Equation Modeling*,
*22*(2), 169-177.

## See also

[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md).

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
bivariate_residuals(fit)
#> =========================================================
#>                BIVARIATE RESIDUALS                       
#> =========================================================
#> Pearson chi-square per item pair, divided by its df.
#> Ranks which pairs strain the model. NOT a calibrated
#> test: referred to chi-square this statistic almost
#> never rejects, so a low value is not evidence of fit.
#> Use n_reps for bootstrap p-values.
#> 
#>          Item1    Item2    Item3    Item4    Item5
#> Item2   0.0422
#> Item3   0.0119   0.0004
#> Item4   0.6329   1.1869   0.4269
#> Item5   0.0702   0.3651   0.0733   0.0279
#> Item6   0.0160   0.2685   0.0699   0.5227   0.1115
#> 
#> Largest: Item4 x Item2 = 1.1869
#> Total BVR: 3.8263
#> =========================================================
```
