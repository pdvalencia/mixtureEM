# Plots for a Fitted Latent Transition Model

Three views, all in base graphics:

- `"prevalence"`:

  latent status prevalence across occasions - the summary of where
  people are over time;

- `"transitions"`:

  a shaded matrix of transition probabilities, one panel per pair of
  adjacent occasions, with the values printed in;

- `"profiles"`:

  item-response probabilities (or means) by status, which is what the
  status labels rest on.

With more than one latent class, `"prevalence"` draws one panel per
class - the classic longitudinal profile plot, in which a stayer class
is the flat one - and `"transitions"` one panel per class per pair of
occasions. `class` restricts either to a single class.

## Usage

``` r
# S3 method for class 'lta_model'
plot(
  x,
  type = c("prevalence", "transitions", "profiles"),
  main = NULL,
  status_labels = NULL,
  colors = NULL,
  class = NULL,
  ...
)
```

## Arguments

- x:

  An object returned by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- type:

  Which view to draw.

- main:

  Plot title.

- status_labels:

  Optional labels for the latent statuses.

- colors:

  Optional colour vector, recycled across statuses.

- class:

  Optional latent class, for a model fitted with `n_classes` \> 1 or
  `mover_stayer = TRUE`.

- ...:

  Ignored.

## Value

`x`, invisibly.
