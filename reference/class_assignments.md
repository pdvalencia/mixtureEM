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
  type = c("modal", "posterior", "both", "viterbi"),
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
  columns; `"viterbi"`, for an `lta_model` only, globally decodes the
  single most probable status sequence (see Details).

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

`type = "viterbi"` decodes globally instead: it returns the single most
probable *sequence* of statuses across all occasions at once, rather
than the most probable status at each occasion taken separately. The two
can disagree — the occasion-by-occasion assignment above may string
together a sequence the model itself gives zero probability, if a
locally-favoured status at one occasion can only be reached by a
transition the model forbids or scores as unlikely, while Viterbi
decoding can never do that. Its `occasion` argument still selects one
occasion's column of that path; it does not change which path is
decoded. With `n_classes` \> 1 the class and the path are decoded
jointly, and the class each case was assigned to is available as
`attr(result, "class_assigned")`.

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
#> 1     1   0.9467880 0.9467880 0.05321203
#> 2     1   0.8100663 0.8100663 0.18993375
#> 3     2   0.8449567 0.1550433 0.84495674
#> 4     2   0.5582606 0.4417394 0.55826060
#> 5     1   0.9030628 0.9030628 0.09693725
#> 6     2   0.7404903 0.2595097 0.74049035
# To relate the classes to an external variable, do not regress on the
# assigned class - use the bias-adjusted third step instead:
# add_outcome(fit, y)
```
