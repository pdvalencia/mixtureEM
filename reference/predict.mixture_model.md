# Classify new cases from a fitted model

Computes the posterior class probabilities of cases that were not in the
estimation sample, holding every parameter of the fitted model fixed.
This is the E-step of the fitted model applied to new rows: by Bayes'
rule a case's posterior is its class-specific likelihood times its prior
class probability, normalised over classes. Vermunt and Magidson (2025)
show that for the standard emissions this is exactly a logistic (or,
with class-specific variances, quadratic) scoring equation, and that
missing items change only its constant; running the E-step gives the
same numbers without assembling the equation, missing values included.

## Usage

``` r
# S3 method for class 'mixture_model'
predict(
  object,
  newdata,
  type = c("modal", "posterior", "both"),
  group = NULL,
  ...
)

# S3 method for class 'lta_model'
predict(
  object,
  newdata,
  type = c("modal", "posterior", "both"),
  occasion = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model from
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md),
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  or
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- newdata:

  A data frame or matrix holding the indicators (and the predictors, for
  a model that has them) of the cases to classify; for a growth model,
  one column per occasion (see Details).

- type:

  What to return, as in
  [`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md):
  `"modal"` the most probable class, `"posterior"` the matrix of
  posterior probabilities, `"both"` a data frame of the two.

- group:

  For a multiple-group model, the new cases' groups: a vector with one
  entry per row of `newdata`, or the name of a column of `newdata`.
  Ignored otherwise.

- ...:

  Unused.

- occasion:

  For a latent transition model, the index of a single occasion to
  return, as in
  [`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md).

## Value

One row (or element) per row of `newdata`, in the format
[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
returns; `NA` for a case that cannot be scored.

## Details

**Matching the new data to the fit.** Indicators are matched to the
fitted model's columns by name, so `newdata` may hold them in any order
and may carry other columns. A factor or character indicator is coded
against the levels it had when the model was fitted, never against the
levels in `newdata`, and a value the fitted model never saw is an error
rather than a silent mis-coding. A binary indicator that was recoded to
0/1 at fit time is recoded the same way. The same holds for a factor or
character predictor: a level the model was not fitted with is an error,
where scoring it would silently treat the case as belonging to the
reference level.

**Growth models.** For a
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
fit, `newdata` holds one row per case and one column per occasion, in
the order of the wide matrix the model was fitted on (time by time, and
within a time item by item for several items). Columns named as the
model's `predictors` are taken as the predictors and every other column
as an occasion, so their number must equal the fitted model's. A model
fitted from long data is scored from the same wide layout. A missing
occasion drops out of the case's likelihood, as in the fit. A
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
fit with `growth_predictors` needs those columns in `newdata` too, by
name (a factor as the factor), since they move each case's own class
means; a case missing one gets `NA`.

**Repeated-measures models.** A
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
fit is scored the same way: one column per item and occasion, in the
order of the wide matrix it was fitted on (time by time, item by item
within a time), plus any predictors by name. The columns are numeric
codes, as the model read them; a binary item the fit recoded to 0/1 (say
from 1/2) is recoded the same way.

**Multiple-group models.** For a model fitted with `group`, each new
case is scored with its own group's parameters: its group's class sizes
when the prevalences differ by group, and its group's item parameters
when the measurement does. Give the groups as `group`, either a vector
with one entry per row of `newdata` or the name of a `newdata` column,
on the values the model was fitted with; a group the model never saw is
an error, and a case with a missing group gets `NA`. A model fitted with
`group_effects = "none"` ignores the group, so `group` is not needed.

**Missing values.** A missing answer simply drops out of that case's
likelihood (full-information scoring), whether or not the estimation
sample had any. A case missing every indicator has no information to
classify on and gets `NA`, the same rule that dropped such rows from the
estimation.

**Predictors and distal outcomes.** For a model with `predictors` (one
step or three steps), the prior class probabilities are the fitted
multinomial regression evaluated at the new case's predictors, so
`newdata` must carry those columns; a factor predictor is supplied as
the factor itself. A case with a missing predictor gets `NA`. A distal
outcome is never used: a new case is classified from its indicators
only, which is what keeps the outcome out of the classification (Vermunt
and Magidson, 2025, sec. 1.4). For a model fitted with a distal outcome
in one step, [`predict()`](https://rdrr.io/r/stats/predict.html) on the
estimation data therefore differs from
[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md),
which conditions on the outcome too. It differs for a three-step model
with predictors as well: there
[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
reports the step-3 posterior, which conditions on the step-1 class
assignment, while [`predict()`](https://rdrr.io/r/stats/predict.html)
conditions on the indicators themselves, since a new case has no step-1
assignment. On a one-step model the two agree exactly.

**Latent transition models.** A
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
fit is scored from one column per item and occasion, as a
repeated-measures model is, and returns a status per occasion in
[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)'s
formats (`occasion` picks one). Each occasion's posterior conditions on
every occasion the case was observed at, earlier and later alike (the
smoothed probability of the forward-backward recursion), which is what
[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
reports for the estimation sample; a mixture of chains, such as a
mover-stayer model, is mixed over its classes the same way. Fits with
predictors, groups, a random intercept, direct effects on items, or
three steps are not covered yet.

**Not covered.** Models with direct effects of covariates on items
(`predictors_items`) and models combining predictors with a distal
outcome in one structural model are refused with an error.

## References

Vermunt, J. K., & Magidson, J. (2025). Linear logistic scoring equations
for latent class and latent profile models: A simple method for
classifying new cases. *Structural Equation Modeling*, *32*(3), 541–549.
[doi:10.1080/10705511.2024.2375736](https://doi.org/10.1080/10705511.2024.2375736)

## See also

[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
for the cases in the estimation sample.

## Examples

``` r
set.seed(1)
X   <- matrix(rbinom(600, 1, rep(c(0.2, 0.8), each = 300)), ncol = 6,
              dimnames = list(NULL, paste0("u", 1:6)))
fit <- fit_mixture(X[1:80, ], n_classes = 2, measurement = "binary")
new <- X[81:100, ]
new[1, 2] <- NA                       # a missing answer is fine
predict(fit, new)
#>  [1] 1 1 2 1 2 1 1 1 2 1 1 1 1 2 1 1 1 2 2 1
head(predict(fit, new, type = "both"))
#>   class probability   Class 1    Class 2
#> 1     1   0.8380779 0.8380779 0.16192215
#> 2     1   0.8842124 0.8842124 0.11578763
#> 3     2   0.8790902 0.1209098 0.87909022
#> 4     1   0.9919797 0.9919797 0.00802027
#> 5     2   0.7016227 0.2983773 0.70162270
#> 6     1   0.9104523 0.9104523 0.08954769
```
