# Class Assignments for Each Case

The class each case is assigned to, and the posterior probabilities
behind that assignment. This is the accessor for the per-case
classification, so that reaching it does not mean indexing into the
fitted object's internals.

## Usage

``` r
# S3 method for class 'lta_model'
class_assignments(
  object,
  type = c("modal", "posterior", "both"),
  occasion = NULL,
  ...
)

class_assignments(object, type = c("modal", "posterior", "both"), ...)

# Default S3 method
class_assignments(object, type = c("modal", "posterior", "both"), ...)
```

## Arguments

- object:

  A fitted model: a `mixture_model` (including the growth models) or an
  `lta_model`.

- type:

  What to return. `"modal"` (default) gives the assigned class;
  `"posterior"` the full matrix of posterior probabilities; `"both"` a
  data frame carrying the assignment, its probability, and the posterior
  columns.

- occasion:

  For an `lta_model`, the index of a single occasion. Omit for every
  occasion, which `type = "both"` does not support.

- ...:

  Passed to methods.

## Value

For `"modal"`, an integer vector of length n. For `"posterior"`, an
n-by-K matrix with the class labels as column names. For `"both"`, a
data frame with `class`, `probability` (the assigned class's posterior
probability, i.e. a per-case classification certainty), and then the K
posterior columns.

## Details

For a latent transition model the assignment is of latent *status*, and
a status assignment is made at every occasion rather than once per case
— the same convention
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md)
and the entropy in
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)'s
`metrics` already follow. Supply `occasion` to work with one occasion at
the shape the mixture methods return; omit it for all of them at once.
To assign the latent *class* of a mixture latent Markov model, use
`object$class_posterior`.

Modal class assignment discards classification error. Do not use the
returned class as though it were an observed variable in a subsequent
regression, ANOVA or t-test: doing so attenuates the association,
severely when the classes are not well separated (Bolck, Croon &
Hagenaars, 2004; Vermunt, 2010; Bakk, Tekle & Vermunt, 2013). Use
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
and
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md),
which correct for it. This function is for plotting, exporting and
describing a solution.

## References

Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
structure models with categorical variables: One-step versus three-step
estimators. *Political Analysis*, *12*(1), 3–27.
[doi:10.1093/pan/mph001](https://doi.org/10.1093/pan/mph001)

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, *18*(4), 450–469.
[doi:10.1093/pan/mpq025](https://doi.org/10.1093/pan/mpq025)

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*,
*43*(1), 272–311.
[doi:10.1177/0081175012470644](https://doi.org/10.1177/0081175012470644)

## See also

[`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md),
[`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md),
[`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md).

## Examples

``` r
set.seed(1)
X   <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
table(class_assignments(fit))
#> 
#>  1  2 
#> 56 44 
head(class_assignments(fit, "both"))
#>   class probability   Class 1    Class 2
#> 1     1   0.9469587 0.9469587 0.05304134
#> 2     1   0.8106031 0.8106031 0.18939689
#> 3     2   0.8447412 0.1552588 0.84474124
#> 4     2   0.5578654 0.4421346 0.55786544
#> 5     1   0.9033559 0.9033559 0.09664411
#> 6     2   0.7401666 0.2598334 0.74016660
# To relate the classes to an external variable, do not regress on the
# assigned class - use the bias-adjusted third step instead:
# add_outcome(fit, y)
```
