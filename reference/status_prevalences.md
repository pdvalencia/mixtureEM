# Latent Status Prevalences by Occasion

The proportion in each latent status at each occasion. The model-implied
prevalences propagate \\\delta\\ through the transition matrices; the
empirical ones average the posterior status probabilities.

With more than one latent class the whole-sample prevalences are
returned by default - the mixture's own marginal, which is what should
be compared with the observed proportions - and `class` picks out a
single class's chain.

## Usage

``` r
status_prevalences(object, type = c("model", "posterior"), class = NULL)
```

## Arguments

- object:

  An object returned by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- type:

  `"model"` (default) or `"posterior"`.

- class:

  Optional latent class, for a model fitted with `n_classes` \> 1 or
  `mover_stayer = TRUE`.

## Value

An occasions-by-statuses matrix.
