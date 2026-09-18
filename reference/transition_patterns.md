# Joint Latent-Status Pattern Table

The joint distribution of latent status across every occasion at once:
one row per possible sequence of statuses, with the count and proportion
of cases following it. This is one of the standard LTA reports and is
not derivable from
[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md)
or
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md),
which each describe only one or two occasions at a time.

## Usage

``` r
transition_patterns(
  object,
  type = c("model", "posterior", "modal"),
  class = NULL
)
```

## Arguments

- object:

  An object returned by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- type:

  `"model"` (default) propagates \\\delta\\ through the transition
  matrices. With `n_classes` \> 1 the classes are summed, weighted by
  `class_weights`, unless `class` selects one. `"posterior"` sums each
  case's actual posterior probability over every possible path; it is
  exact but there is no shortcut for it – the path posterior does not
  factorise, so every path is enumerated directly – and it is
  additionally refused for a model with more than one latent class,
  where `"modal"` is the alternative. `"modal"` cross-tabulates the
  joint-MAP (Viterbi) path from `class_assignments(object, "viterbi")`.
  That is a different thing from a cross-tabulation of the per-occasion
  modal statuses, which can put mass on a pattern the model itself gives
  zero probability; see
  [`class_assignments()`](https://pdvalencia.github.io/mixtureEM/reference/class_assignments.md)'s
  own documentation of the distinction. A per-occasion modal
  classification table is the more common summary of a latent transition
  model, not this joint decode – get it with
  `table(class_assignments(object, "modal"))` rather than
  `type = "modal"` if that is the number being matched against.

- class:

  Optional latent class, for a model fitted with `n_classes` \> 1.
  Applies only to `type = "model"`; ignored otherwise.

## Value

A data frame with one integer column per occasion, named from the
model's time labels, then `count` and `proportion`.

## Details

The table has one row per possible status sequence, so it has \\K^T\\
rows whatever `type` is asked for, and every `type` is refused above
\\K^T \> 10000\\. No `type` escapes that: the enumeration is the shape
of the answer, not a way of computing it. To describe a model with more
occasions than that allows, use
[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md)
or
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md),
which look at one or two occasions at a time, or decode each case's own
path with `class_assignments(object, "viterbi")`, which returns one row
per *case* and so is bounded by the sample rather than by \\K^T\\.

## See also

[`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md),
[`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md).
