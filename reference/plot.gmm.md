# Trajectory Plot for a Growth Mixture Model

Draws the estimated mean trajectory of each class, with the observed
data of the cases assigned to it behind the curves.

## Usage

``` r
# S3 method for class 'gmm'
plot(
  x,
  observed = c("means", "cases", "none"),
  main = NULL,
  class_labels = NULL,
  colors = NULL,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md).

- observed:

  What of the observed data to draw: `"means"` (the observed mean at
  each occasion among cases modally assigned to the class, dotted),
  `"cases"` (individual trajectories, translucent) or `"none"`.

- main:

  Plot title.

- class_labels:

  Optional class labels for the legend.

- colors:

  Optional colour vector, recycled across classes.

- ...:

  Ignored.

## Value

`x`, invisibly.
