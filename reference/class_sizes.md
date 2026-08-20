# Class Sizes of a Fitted Mixture Model

Returns the estimated size of each latent class in the three forms
applied papers report: the model's class proportion, the expected number
of cases, and the number of cases modally assigned to the class. Case
weights are used when the model was fitted with any.

## Usage

``` r
class_sizes(object, ...)

# S3 method for class 'mixture_model'
class_sizes(object, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Passed to methods.

## Value

A data frame with one row per class: `class`, `proportion`
(model-estimated class weight), `n_expected` (proportion times the
analysed sample size), and `n_modal` (cases assigned by highest
posterior probability).

## Details

The smallest class is one of the things readers judge a solution by, and
there are two published conventions for it. Lee et al. (2023, p. 654):
"If the smallest class contains less than 5\\ size for the smallest
class is less than 25, it is recommended that the model only be retained
as the optimal model if the researcher can accurately defend what is
gained from this small class given the possibility of low power and a
lack of statistical precision." Jung and Wickrama (2008, p. 312) give a
weaker floor among their checks: no less than 1\\

Both are reporting conventions, and mixtureEM enforces neither. This
function applies no threshold, raises no warning and filters nothing; a
small class is estimated and returned like any other. What the
conventions ask of you is a defence, not a deletion: keep the class if
it can be justified substantively, drop it if it cannot, and report the
number either way.

## References

Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class
growth analysis and growth mixture modeling. *Social and Personality
Psychology Compass*, *2*(1), 302-317.
[doi:10.1111/j.1751-9004.2007.00054.x](https://doi.org/10.1111/j.1751-9004.2007.00054.x)

Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction
to growth mixture models (GMM). In *International Encyclopedia of
Education* (4th ed., Vol. 14, pp. 646-655). Elsevier.
[doi:10.1016/B978-0-12-818630-5.10076-4](https://doi.org/10.1016/B978-0-12-818630-5.10076-4)

## See also

[`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)
for the per-case assignment the `n_modal` column counts.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
class_sizes(fit)
#>   class proportion n_expected n_modal
#> 1     1  0.5844398   58.44398      56
#> 2     2  0.4155602   41.55602      44
```
