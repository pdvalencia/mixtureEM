# Print Classification Diagnostics

Prints the two tables that describe how cleanly the model assigns cases.

The **Average Posterior Probability (AvePP)** matrix has one row per set
of observations modally assigned to a class and one column per class,
holding the mean posterior probability. High values on the diagonal, low
values off it, indicate well-separated classes.

The **classification table** cross-classifies the probabilistic
memberships against the modal assignment, and yields the classification
error: the proportion of cases the modal rule is expected to place in
the wrong class. See
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md)
for the details, and
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
and
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
for the fit of the model itself rather than the quality of its
assignments.

The diagonal of the AvePP matrix reads as a percentage of the cases
assigned to that class: a value of 0.91 "suggests that 91% of subjects
in the assigned class fit that category, while 9% of the subjects in
that class do not accurately fit that category" (Fanti & Henrich, 2010,
as reported by Lee et al., 2023, p. 653).

This function does not compute entropy, which is the other number
solutions are judged by. That is `fit$metrics$entropy`, printed as
`Rel. Entropy` and carried as the `Entropy` column of
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md),
where it is documented.

Both tables use the case weights when the model was fitted with any.

Three further blocks print when they apply. None carries a threshold, a
warning or a verdict; each is a measurement to report and judge.

**Class separation** (all-continuous models only): the Mahalanobis
distance between each pair of class means, against the pooled
within-class variance of each indicator. It is a property of the
parameters on the indicators' own scale, not of how cleanly cases were
assigned, which is what entropy measures; with one indicator it is
Cohen's d. The smallest distance is the one that limits class
enumeration: Tein, Coxe and Cham (2013) found "little power to detect
the correct number of classes for small (d = .2) or even medium (d = .5)
inter-class distances, under any conditions or with any selection
methods", and that sample size mattered far less than separation. Read
it in one direction only. When the true classes are poorly separated the
fitted classes are not the true ones, and their distance overstates the
truth: in a simulation of three classes at true neighbouring distance 1
(n = 600) the reported smallest distance averaged 1.6, while at true
distance 3 it averaged 2.9. A small value says the classes are close; a
moderate one does not by itself show they are not. For categorical
indicators no distance is reported, since Cohen's d does not exist
there.

**Modal minus probability-weighted class means** (continuous and binary
indicators; a binary mean is the proportion endorsing): each class's
indicator means computed over the cases modally assigned to it, minus
the same means weighted by the posterior probabilities of the
measurement model. It shows how much exporting a modal class variable
and analysing it as if it were known changes what the classes look like
(Chen et al., 2026). The package's own answer to classification error is
to correct for it, with
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
and
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md);
this block is for the analysis that will use the exported variable
anyway.

**Cases moved by the third step** (a structural model attached with
`correction = "ML"`): the modal class under the measurement model's
posterior against the modal class under the step-3 posterior, as a
count, a proportion and a cross-tabulation. The three-step method
assumes the structural model leaves the classes as step 1 defined them
(Vermunt, 2010; Asparouhov & Muthén, 2014); many moved cases mean the
structural model is redefining the classes, and the one-step and
two-step fits are the comparison to run.

## Usage

``` r
classification_diagnostics(object, ...)

# Default S3 method
classification_diagnostics(
  object,
  n_boot = 0,
  level = 0.95,
  random_state = 123,
  n_cores = .default_n_cores(),
  ...
)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Passed to methods.

- n_boot:

  Non-negative integer. Number of case-resampling bootstrap draws behind
  a percentile interval for each class proportion; `0` (the default)
  skips it. Each draw continues the fitted solution on the resampled
  cases (one EM run from the fitted parameters, no restarts), so its
  classes keep the fitted labels. Available for an unconditional,
  unweighted, single-group, one-step fit. The interval can be too narrow
  when some item probabilities sit at 0 or 1: in a simulation from such
  a four-class solution it covered the true proportion only 69\\ the
  time for the classes involved.

- level:

  Confidence level of that interval. Default `0.95`.

- random_state:

  Integer seed for the resampling. Default `123`.

- n_cores:

  Positive integer. Processes to spread the draws over. Default `1`, or
  `options(mixtureEM.n_cores = )` where set.

## Value

Invisibly, a list with `ave_pp` (the K x K matrix), `table` (the
classification table), `error` (the classification error) and `classes`,
Masyn's (2013) per-class table: `Proportion` (the model's class
proportion), `Lower` and `Upper` (its bootstrap interval, when
`n_boot > 0`), `mcaP` (the share of cases modally assigned to the
class), `AvePP` (the diagonal of the AvePP matrix) and `OCC`, the odds
of correct classification, \\\[\mathrm{AvePP}\_k/(1-\mathrm{AvePP}\_k)\]
/ \[\pi_k/(1-\pi_k)\]\\: how much better modal assignment does than
guessing from the class proportion alone. Masyn suggests values above 5
indicate good separation. Also `separation`, the K x K distance matrix
with its smallest off-diagonal value as attribute `"min"`;
`profile_shift`, the K x indicators matrix of modal minus
probability-weighted means; and `class_shift`, a list of `table`,
`n_moved` and `prop_moved`. Each is `NULL` where it does not apply. All
are also printed to the console.

## References

Asparouhov, T., & Muthén, B. (2014). Auxiliary variables in mixture
modeling: Three-step approaches using Mplus. *Structural Equation
Modeling*, *21*(3), 329-341.
[doi:10.1080/10705511.2014.915181](https://doi.org/10.1080/10705511.2014.915181)

Celeux, G., & Soromenho, G. (1996). An entropy criterion for assessing
the number of clusters in a mixture model. *Journal of Classification*,
*13*(2), 195-212.
[doi:10.1007/BF01246098](https://doi.org/10.1007/BF01246098)

Chen, X., Chen, S., Wang, C., Sweeney, C., Bailey, R., & Samsudin, N.
(2026). When is hard class assignment defensible? An uncertainty-aware
framework for psychometric profile interpretation. *Educational and
Psychological Measurement*. Advance online publication.
[doi:10.1177/00131644261469071](https://doi.org/10.1177/00131644261469071)

Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction
to growth mixture models (GMM). In *International Encyclopedia of
Education* (4th ed., Vol. 14, pp. 646-655). Elsevier.
[doi:10.1016/B978-0-12-818630-5.10076-4](https://doi.org/10.1016/B978-0-12-818630-5.10076-4)

Nagin, D. S. (2005). *Group-Based Modeling of Development*. Harvard
University Press.

Tein, J.-Y., Coxe, S., & Cham, H. (2013). Statistical power to detect
the correct number of classes in latent profile analysis. *Structural
Equation Modeling*, *20*(4), 640-657.
[doi:10.1080/10705511.2013.824781](https://doi.org/10.1080/10705511.2013.824781)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450-469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

Masyn, K. E. (2013). Latent class analysis and finite mixture modeling.
In T. D. Little (Ed.), *The Oxford Handbook of Quantitative Methods*
(Vol. 2, pp. 551-611). Oxford University Press.

## See also

[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
for the per-case assignment these diagnostics summarise.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
classification_diagnostics(fit)
#> =========================================================
#>           AVERAGE POSTERIOR PROBABILITIES (AvePP)        
#> =========================================================
#> Rows: Modal Assignment | Columns: Mean Probability
#> 
#>                  Prob C 1 Prob C 2
#> Assigned Class 1    0.823    0.177
#> Assigned Class 2    0.282    0.718
#> =========================================================
#> 
#> =========================================================
#>                CLASSIFICATION TABLE                      
#> =========================================================
#> Rows: model-expected membership | Columns: modal assignment
#> 
#>         Modal 1 Modal 2    Total
#> Class 1 46.0756 12.4153  58.4908
#> Class 2  9.9244 31.5847  41.5092
#> Total   56.0000 44.0000 100.0000
#> 
#> Classification error: 0.2234 (22.34% of 100 cases)
#> =========================================================
#> 
#> Class sizes and classification quality (Masyn, 2013)
#>  Class Proportion mcaP AvePP   OCC
#>      1      0.584 0.56 0.823 3.306
#>      2      0.416 0.44 0.718 3.572
#> 
#> Modal minus probability-weighted class means of the indicators 
#> (how far exporting a modal class variable moves each profile)
#>             V1     V2     V3     V4     V5
#> Class 1 -0.179 -0.021  0.016  0.023  0.019
#> Class 2  0.197  0.031 -0.036 -0.020 -0.018
```
