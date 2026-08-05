# Profile Plot for a Fitted Mixture Model

Draws a profile plot of the measurement model: one line per latent
class, with every indicator placed on a common \[0, 1\] axis. Binary
indicators are shown as endorsement probabilities; continuous indicators
are min-max scaled against their observed range; polytomous indicators
are summarised by their expected category and scaled to \[0, 1\].
Rescaled items are marked with "\*".

Uses only base graphics and the colour-blind-friendly Okabe-Ito palette.

## Usage

``` r
# S3 method for class 'mixture_model'
plot(
  x,
  main = "Latent Class / Profile Plot",
  class_labels = NULL,
  colors = NULL,
  ...
)
```

## Arguments

- x:

  A fitted `mixture_model` object.

- main:

  Plot title.

- class_labels:

  Optional character vector of labels for the classes. Defaults to
  "Class 1", "Class 2", ...

- colors:

  Optional vector of colours (one per class). Defaults to the Okabe-Ito
  palette, recycled if necessary.

- ...:

  Currently unused; present for S3 compatibility.

## Value

The fitted model, invisibly.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
plot(fit)

```
