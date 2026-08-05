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

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2)
class_sizes(fit)
#>   class proportion n_expected n_modal
#> 1     1  0.5835382   58.35382      56
#> 2     2  0.4164618   41.64618      44
```
