# Trajectory Plot for a Latent Class Growth Model

Draws the estimated trajectory of each class on the response scale, with
the occasions on the x-axis — the figure an LCGA is read from, since the
classes *are* the trajectories — and the observed data of the cases
assigned to each class behind them.

## Usage

``` r
# S3 method for class 'lcga'
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
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md).

- observed:

  What of the observed data to draw: `"means"` (the observed mean at
  each occasion among cases modally assigned to the class, dotted),
  `"cases"` (individual trajectories, translucent — informative for a
  continuous or count outcome, much less so for a binary one) or
  `"none"`.

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
