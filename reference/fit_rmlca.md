# Repeated-Measures Latent Class Analysis

Fits a repeated-measures latent class model to data in which the same
indicators are observed at several occasions. Each person belongs to one
latent class for the whole study, so the classes describe
*trajectories*: patterns of response that span the occasions rather than
a snapshot at one of them. This is Collins and Lanza's RMLCA (2010, sec.
7.2); with continuous indicators the same model is Wang and Wang's
longitudinal latent profile analysis (2020, sec. 6.3.1), and it is
obtained here simply by setting `measurement = "continuous"`.

The likelihood is that of an ordinary latent class model applied to the
\\J \times T\\ stacked indicators, \$\$P(y_i) = \sum_k \gamma_k \prod_t
\prod_j \rho\_{jtk}(y\_{ijt}),\$\$ so everything the package already
offers — model selection, the bootstrap likelihood-ratio test,
predictors of class membership, distal outcomes, survey designs and FIML
for missing data — applies unchanged.

`measurement_invariance` controls whether the item-response parameters
\\\rho\_{jtk}\\ are held equal across occasions. Constraining them makes
a class label mean the same thing at every occasion and sharply reduces
the number of free parameters; leaving them free lets an item behave
differently over time. The two models are nested, so
[`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
tests the restriction directly.

## Usage

``` r
fit_rmlca(
  indicators,
  n_classes = 2,
  times = NULL,
  measurement = "binary",
  measurement_invariance = c("none", "full", "partial"),
  invariant_items = NULL,
  layout = c("time_major", "item_major"),
  id = NULL,
  time = NULL,
  items = NULL,
  item_names = NULL,
  time_labels = NULL,
  predictors = NULL,
  ...
)
```

## Arguments

- indicators:

  The repeated indicators. Either a wide matrix or data frame with \\J
  \times T\\ columns (see `layout`), a three-dimensional array with
  dimensions n by items by times, or a long data frame together with
  `id` and `time`.

- n_classes:

  Integer. Number of latent classes.

- times:

  Integer. Number of occasions. Required for wide input; inferred
  otherwise.

- measurement:

  Measurement model for a single occasion's items: `"binary"`,
  `"categorical"`, `"continuous"`, or a named list for a mixed block (as
  in
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)).

- measurement_invariance:

  Whether the item parameters are held equal across occasions. `"none"`
  (the default) estimates them separately at every occasion, which is
  usually what you want here: the classes are patterns of change, so
  forcing the items to behave identically over time can erase the very
  differences being modelled. `"full"` holds every item equal, and
  `"partial"` holds only the items named in `invariant_items`.

- invariant_items:

  Item indices or names held equal across occasions. Used only when
  `measurement_invariance = "partial"`.

- layout:

  For wide input, whether columns run `"time_major"` (all items of
  occasion 1, then all items of occasion 2, ...) or `"item_major"` (all
  occasions of item 1, then all occasions of item 2, ...).

- id, time:

  For long input, the case and occasion identifiers, given either as
  column names or as vectors.

- items:

  For long input, the columns to treat as indicators.

- item_names, time_labels:

  Optional display labels.

- predictors:

  Optional predictors of class membership. A grouping variable (Collins
  and Lanza, sec. 7.2.1) is entered this way: a multiple-group model and
  a model with the group as a dummy predictor are equivalent when
  measurement is invariant across groups (their sec. 6.10.2).

- ...:

  Further arguments passed to
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  such as `outcome`, `n_init`, `random_state`, `weights`, `strata` or
  `cluster`.

## Value

An object of class `c("rmlca", "mixture_model")`. In addition to the
usual fields it carries `$longitudinal`, holding the item and occasion
labels, the invariance specification and the wave-missingness pattern.

## See also

[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
for a model in which class membership may change between occasions, and
[`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
for testing invariance.
