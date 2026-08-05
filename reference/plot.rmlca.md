# Trajectory Plot for a Repeated-Measures Latent Class Model

Draws one panel per item with occasions on the x-axis and one line per
latent class, which is the natural way to read RMLCA output: the classes
*are* the trajectories. Pass `type = "profile"` for the single-panel
indicator profile used by
[`plot.mixture_model()`](https://pdvalencia.github.io/mixtureEM/reference/plot.mixture_model.md).

## Usage

``` r
# S3 method for class 'rmlca'
plot(
  x,
  type = c("trajectory", "profile"),
  main = NULL,
  class_labels = NULL,
  colors = NULL,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).

- type:

  `"trajectory"` (default) or `"profile"`.

- main:

  Plot title.

- class_labels:

  Optional class labels for the legend.

- colors:

  Optional colour vector, recycled across classes.

- ...:

  Passed to the profile plot when `type = "profile"`.

## Value

`x`, invisibly.
