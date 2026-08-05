# Transition Probability Matrices

Returns the estimated transition probabilities \\\tau\_{s_t \|
s\_{t-1}}\\: the probability of occupying each latent status at one
occasion given the status held at the previous one. Rows are the origin
status and sum to one; the diagonal is the probability of staying put.

There is one matrix per pair of adjacent occasions unless the model was
fitted with `tau_homogeneous = TRUE`, in which case a single matrix is
shared.

In a mixture latent Markov model each latent class has its own
transitions, so the result gains an outer level indexed by class;
`class` picks one out.

## Usage

``` r
transition_matrix(object, occasion = NULL, class = NULL)
```

## Arguments

- object:

  An object returned by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- occasion:

  Optional index of a single transition (1 for the move from occasion 1
  to occasion 2, and so on). Omit to get them all.

- class:

  Optional latent class, for a model fitted with `n_classes` \> 1 or
  `mover_stayer = TRUE`. Omit to get every class.

## Value

A matrix, or a named list of matrices when `occasion` is omitted and the
transitions are time-heterogeneous, nested inside a list over classes
when the model has more than one.
