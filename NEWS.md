# mixtureEM (development version)

## `measurement_summary()` now shows the sample marginal

Every table gains an `Overall` column between the indicator name and the
classes, and the returned data frame gains an `overall` column beside
`estimate`. It is the observed marginal for that item: the weighted sample
proportion beside a probability, the weighted sample mean beside a mean or a
rate, and for a polytomous item the share of cases in each category.

A conditional number is not readable on its own. A class endorsing an item at
.62 is unremarkable where the sample sits at .60 and is most of what defines the
class where the sample sits at .12, and the table could not tell those apart —
the reader had to go back to the data. Putting the marginal first in the row
makes each class parameter read as a departure from it.

The benchmark is the *observed* marginal rather than the model-implied one,
because a model-implied column would agree with the class parameters by
construction and so could not act as a check on them. It uses the case weights
where the fit has any, and reads the indicators as the fit stored them, which is
after any binary recode — so the proportion is of the same level the probability
beside it is of. For a latent transition model the marginal is the one for the
occasion being printed, or, where the measurement model is held equal across
occasions, the one pooled over all of them.

The column is dropped, with a note saying why, for a fit that does not store its
raw indicators and for a multiple-group measurement model, whose per-group item
parameters cannot be matched to the stacked indicator columns by name; the data
frame then carries `NA` there. Growth models always carry `NA`, a growth factor
mean having no sample marginal to be compared against. No fitted number changes.

## The standardized profile can now be drawn as lines

`plot(fit, type = "line")` draws the same z-scored conditional means as
`type = "bar"`, but as one connected line per class with a zero reference line.
Like the bar chart it needs an all-continuous measurement model, and it takes
the same `scale` argument.

The two are different readings of the same numbers. Bars group by indicator, so
the eye compares classes one indicator at a time; a line follows a single class
across all of them, which is what shows whether two classes differ in level or
in pattern. That is the reading "profile" names in the applied literature, and
until now the only line plot on offer was the default `type = "profile"`, whose
min-max axis has no meaningful origin and whose shape moves with the sample's
most extreme observation. For an all-continuous model, prefer `"line"`.

One helper computes the heights for both renderers, so the two can never
disagree about a number.

## Bivariate residuals can now be calibrated by bootstrap

`bivariate_residuals()` takes `n_reps`, which replaces the statistic's
chi-square reference with a parametric bootstrap and attaches a matrix of
p-values that the print method shows beside each residual.

The reference distribution is the problem being solved. In the simulation of
Oberski, van Kollenburg and Vermunt (2013), a bivariate residual referred to
chi-square rejected at a nominal five percent in zero of two hundred samples in
*every one* of eight null conditions; its empirical mean was between 0.25 and
0.36 where the reference has 1. A statistic that never rejects under a true
model also says little under a false one, so a low bivariate residual is not
evidence of good fit — and the ranking it supports is the least powerful of the
three methods those authors compared. The documentation now says all of this
plainly, along with two limits it had not stated: power falls as the classes
separate, and under a missing-at-random mechanism a large residual is ambiguous
between local dependence and selection bias.

The default `n_reps = 0` is the previous behaviour at the previous cost, and
returns bit-identical residuals. `100` is the recommended working value and
`500` the publication-grade one.

## A standardized profile bar chart

`plot(fit, type = "bar")` draws the grouped bar chart applied latent-profile
papers publish: indicators along the x-axis, one bar per class within each
group, and a zero line separating above-average from below-average conditional
means. It requires an all-continuous measurement model; `type` defaults to
`"profile"`, so the existing figure is unchanged.

Bar heights are z-scores rather than the profile plot's min-max scaling, which
is hostage to a single extreme observation and has no meaningful origin. A
z-score is scale-free, so the figure comes out the same whether or not the
indicators were standardized before fitting — standardizing becomes a display
choice rather than a step in preparing the data. `scale = "within"` divides by
the model-implied within-class standard deviation instead, giving a
Cohen's-d-like reading against residual rather than total dispersion.

## Step-3 standard errors: which one to compare against another program

The `se` documentation now records that `"corrected"`, the default, is the
statistically right answer, while `"robust"` is the *comparability* setting.
Another program reports the step-3 sandwich alone, so reproducing its standard
errors requires asking for `"robust"`; under the default a user checking
mixtureEM against it sees wider intervals, and that difference is a difference
in estimator — the corrected form carries step-1 uncertainty the sandwich omits
— rather than a bug in either program. No estimates or standard errors change.

## `measurement` is now required

`fit_mixture()`, `compare_mixtures()` and `blrt()` no longer default
`measurement` to `"binary"`. Omitting it is an error that lists the valid
types, shows the mixed-type syntax, and suggests a type read off your columns:

```
`measurement` must be specified. Valid types: "binary", "categorical",
"continuous", "count".
Your 8 indicator columns all take two values, so you probably want
  measurement = "binary"
For mixed types:
  measurement = list(binary = 1:5, continuous = 6:8)
```

The suggestion is a hint to confirm, not a choice the package makes. The
storage mode of a column does not determine its measurement model: a 1-5
column is a legitimate `"categorical"`, `"continuous"` or `"count"` indicator,
and the class solution differs across the three. Inferring the type would
settle a modelling question by inspecting storage mode and would make a
script's meaning depend on the data it is run against; a constant default is
the same guess with the data-dependence removed.

This breaks any call that relied on the default. The fix is to add
`measurement = "binary"`, which reproduces the previous behaviour exactly.
`blrt(from_fit = )` is unaffected, since it reads the specification off the
fitted model.

## Continuous indicators default to equal variances across classes

`fit_mixture(measurement = "continuous")` and `compare_mixtures(measurement =
"continuous")` now fit the homoscedastic latent profile model, holding each
item's variance equal across the classes. This changes the estimates, the fit
indices and the class solution of any continuous fit that did not set
`variances_equal`. `variances_equal = FALSE` recovers the previous behaviour
exactly.

The reason is that the unrestricted normal-mixture likelihood is unbounded:
send a class mean to any single data point and that class's variance to zero
and the likelihood diverges, so no maximum likelihood estimate exists and what
the EM reports is a local optimum. Holding the variances equal bounds the
likelihood, and the constrained estimator is consistent (Day, 1969; Hathaway,
1985). Freeing them also invites classes that describe non-normality in a
single population rather than distinct subgroups (Bauer & Curran, 2003).

The restriction is substantive and testable, and the expectation is that you
fit both and compare. It is not the safe choice but the well-posed one: the
homoscedastic model fails visibly, by splitting a genuinely heteroscedastic
class in two, while the free model fails silently, as a boundary solution that
gets written up as a finding.

Only the two user-facing entry points resolve this default. The growth,
time-block and group-block paths — LCGA, GMM and RMLCA — are unchanged.

## `add_covariates()` and `add_outcome()` accept a formula and `data`

Both functions gain a `data` argument, so the columns can be named instead of
extracted:

```r
add_covariates(fit, ~ T1age + T1sex + T1SHexp, data = df)
add_outcome(fit, ~ T3NHR_1, data = df)
```

`predictors` and `covariates` also accept a character vector of column names.
`add_outcome()`'s formula must name exactly one column, since one call fits one
distal outcome.

This is a matter of typing, and nothing more. The existing calling style —
`add_covariates(fit, df[, c("T1age", "T1sex")])`, `add_outcome(fit, df$T3NHR_1)`
— is unaffected and unchanged, with no deprecation and no message: passing a
computed vector such as `scale(y)` is often the right thing to do, because the
variable is not a column of anything. Both forms meet at the same code as soon
as the columns are in hand, and every estimate is identical either way.

## Two defaults change: the class-membership prior now reaches step 3

`bayes_constants$latent` is documented as the Dirichlet prior on the class
probabilities. It applied when those probabilities were estimated as K-1 free
weights, and silently stopped applying the moment covariates entered and they
became a regression instead. That leak is now closed: the prior applies to the
class-membership regression too, written as fractional pseudo-data — one row per
class per unique covariate pattern, weight `latent / (K * U)`, adding
`latent / K` cases to each class.

Two consequences, both deliberate:

* **Covariate coefficients move.** Every one of them shrinks slightly towards
  zero, because the prior makes the class sizes slightly more equal.
* **Their standard errors move.** They shrink, because unlike the "ghost"
  observation the prior replaces, these rows are part of the objective being
  maximised and so enter the information matrix.

`bayes_constants = list(latent = 0)` restores the previous behaviour exactly, in
both the coefficients and the standard errors; the ghost observation that guards
against complete separation is kept in that case, as before.

## The three-step corrections take an `assignment` argument

`add_covariates()`, `add_outcome()` and `fit_mixture()` gain
`assignment = c("proportional", "modal")`: how step 1's posteriors are turned
into the assigned-class variable whose classification error the BCH and ML
corrections invert. **The default does not change** — `"proportional"` follows
Bakk, Tekle and Vermunt (2013), who found it at least as accurate as modal
assignment across 54 simulation conditions and clearly better when the classes
are poorly separated. Use `assignment = "modal"` when reproducing an analysis
whose classes were assigned that way. The rule in force is stored on the fit and
printed next to the correction, so a saved model still says which one produced
it.

`add_covariates()` also now documents two things that matter when comparing
coefficients with a published set: a case missing a predictor is retained and
completed under the class-invariant Gaussian marginal rather than listwise
deleted, so the analysed N can differ; and `se = "corrected"` carries the step-1
uncertainty where `se = "robust"` reports only the step-3 sampling variability.

## New: `class_assignments()`

The per-case classification now has an accessor, so reaching it no longer means
writing `max.col(fit$log_resp)` by hand. `type = "modal"` gives the assigned
class, `"posterior"` the full matrix, and `"both"` a data frame carrying the
assignment together with its probability — a per-case classification certainty.
It works on `mixture_model`, the growth models, and `lta_model`, where the
assignment is of latent status and an `occasion` argument picks one out.

Its documentation carries the warning that is the reason it exists: do not use
the returned class as though it were an observed variable in a subsequent
regression, ANOVA or t-test. Use `add_covariates()` and `add_outcome()`, which
correct for the classification error that discards.

## `print()` now shows the full set of fit indices

`print()` on a fitted model showed the log-likelihood and relative entropy;
reading its BIC meant reaching into `fit$metrics` or running a one-model
`compare_mixtures()`. It now prints `Log-Likelihood`, `Parameters`, `AIC`,
`BIC`, `SABIC` and `Rel. Entropy` — exactly the columns `compare_mixtures()`
tabulates, so printing one model and comparing a range of K can never show two
different sets of numbers for the same fit. `print()` on a latent transition
model gains the two indices it was missing.

On a three-step fit the criteria are read off the same set of metrics as the
log-likelihood above them, never a mixture of the two. On a fit whose variances
collapsed, the BIC line says so where it appears, since that number is inflated
by the spike and is not comparable with a clean fit's.

## The `n_init` advice now scales with the search that actually ran

The "refit with `n_init = 100`" advice fired whenever the maximum was found by a
single start, whatever `n_init` had been — so a user who ran `n_init = 200` was
told to refit with 100. The advice now reads both of the counts the fit carries,
the restarts requested and the restarts run out to convergence, and says one of
three things.

Below 100 requested it is unchanged. Above 100 requested but with fewer than 100
run out to convergence — which is the ordinary case on a staged search, where
only the most promising restarts are refined — it says that the maximum failed
to replicate among the restarts that were run out, that this is a thinner test
than the requested count makes it sound, and that `n_init` should be raised
further before anything is read into it. Only when 100 or more restarts reached
convergence does it raise the specification: how well separated the classes are,
how heavily parameterised the within-class structure is, and whether there are
more classes than the data support.

The strongest reading is now hedged where it was asserted. The number of random
starts a mixture needs grows with the number of classes, the number of free
parameters and how poorly the classes separate, so the message no longer claims
that more starts are unlikely to help, and no longer ranks over-extraction ahead
of the other causes. This applies to `print()`, the warning, `compare_mixtures()`
and `compare_longitudinal()`, which share one helper. No number changes.

## Diagnostics now say what to change, and by how much

No default changed anywhere in this group, so every existing fit returns the
same numbers. What changed is what the package tells you about them.

* **An unreplicated maximum is now a warning, not a line in `print()`.** The
  most informative local-maximum signal there is — "the best solution was found
  by 1 of 20 starts" — was a `cat()` line, invisible to anyone working from
  `summary()` or from the coefficients. It is now a warning that names the
  remedy: refit with `n_init = 100`, and if the maximum still does not
  replicate, read that as a problem with the specification rather than with the
  search. It fires only from ten requested starts upward, below which a lone
  replication carries no information, and it stays silent on a fit that has
  already been flagged for a collapsed variance or a growth-factor boundary —
  those warnings say that raising `n_init` can make matters worse, and two
  warnings must not give opposite advice about the same argument.

* **The number of restarts is now reported honestly on the staged searches.** A
  growth mixture model at `psi = "equal"`, or a latent transition model with
  more than one class, ranks its restarts on a short pass and carries only three
  of them to convergence. The report counted the survivors, so a 50-start search
  announced itself as a 10-start one; `fit_lta()` kept no counts at all. Both
  numbers are now carried, and printed as "found by 1 of 3 starts that ran to
  convergence (of 50 requested)".

* **`compare_mixtures()` and `compare_longitudinal()` gain an `Unreplicated`
  column**, with a line after the table naming the class counts to refit before
  reporting. The per-model warning is suppressed inside
  `compare_longitudinal()`, which would otherwise raise it once per K before the
  table it belongs next to had been printed.

* **The non-convergence warning names a new `max_iter`** — double the one that
  failed — instead of saying "a larger `max_iter`", and `fit_lta()` now issues
  it at all. It has its own EM driver, so it reported non-convergence only
  through `print()`.

* **`blrt()` warns when 100 draws cannot resolve the decision.** A bootstrap
  p-value can only take the values 1/(B+1), 2/(B+1), …, so at 100 draws it
  cannot separate .04 from .06; when the result lands within one step of .05 the
  test now says so and names `n_reps = 999`. It also counts draws where the
  larger model fitted worse than the model nested inside it — a symptom of a
  replicate search that stopped short — reports the count as `n_negative`, and
  recommends `n_init_boot = 50`.

* **The boundary-probability note now offers a way forward**, rather than
  stating the problem and stopping: read it substantively, since an item every
  member of a class answers identically is often the finding, or refit with a
  stronger `bayes_constants = list(categorical = ...)` if that parameter needs
  an interpretable standard error.

* **Every default the estimator chose is now traceable to a source.**
  `vignette("estimation")` gains three sections and a rewritten one: why
  `n_init = 20` is a floor and what the published replication rates actually
  say about a count of 1 (including the correction that a 3–10% band puts 1 of
  20 *inside* it, which is an argument for more starts rather than a verdict);
  how long EM runs and why the doubling escalation; the bootstrap test, its
  p-value formula, and what 100 draws can and cannot resolve; and a table of
  every warning the package raises with the action for each. Where a number
  cannot be sourced — the `1e-4` EM tolerance, `n_init_boot = 10` — the
  vignette says so rather than dressing it up.

* **The growth-model help files now carry the applied reporting conventions,
  with their sources, and no behaviour changed.** `?fit_gmm` names the three
  specification levels applied papers use, says which of them `psi` and
  `residual_equal` correspond to and which this package cannot fit, and states
  the cost of the `psi = "equal"` default — the same constraint that buys
  stability can buy an extra class that is an artefact of it. `?fit_lcga` says
  why an LCGA's information criteria can improve monotonically with K, and why
  that is a symptom rather than a result. Both now say how many occasions each
  polynomial degree needs. `?class_sizes` gives the two published small-class
  conventions and states plainly that the package enforces neither, and the two
  comparison functions document how to read the `Entropy` column — anchors, and
  the fact that entropy is not evidence for the number of classes.

* **The `categorical` and `latent` priors now carry their evidence too.** The
  documentation justified the `variances` prior and left the other two as bare
  defaults. Both the strength of one added observation and the decision to
  spread it in agreement with each item's observed marginal — rather than
  uniformly over the cells — come from Galindo Garre and Vermunt (2006), whose
  simulation finds that form the most accurate of those studied and shows why a
  uniform spread degrades as the number of items grows. No value changed; the
  package already implemented the prior their results favour.

* **`fit_lta()`'s `bayes_constants` now reaches the measurement model.** The two
  priors this model uses are documented as a division of labour — `smoothing`
  for the status prevalences and the transition matrices, `bayes_constants` for
  the measurement model — and the code did not implement it. `smoothing` was
  passed into the measurement model's M-step as well, where it took precedence,
  so `bayes_constants = list(categorical = ...)` had no effect at all and
  `smoothing = 0` quietly removed the measurement prior too, returning the
  item-response probabilities of exactly 0 and 1 that prior exists to prevent.
  Fits using the defaults are unchanged to every digit, both arguments being 1.

* **`fit_lta()` now reports how much of each transition row comes from the prior
  rather than from the data**, and says so when it exceeds five percentage
  points on any row, naming the row and the size of the effect. An origin status
  that few cases occupy is a row the one pseudo-case of smoothing carries a
  visible share of, and the share has an exact form rather than needing to be
  estimated. It is a reading caution, not a verdict on the fit: those
  transitions should be reported as indicative. The remedy it points at first is
  `transition_invariance = "full"`, which puts every occasion's cases behind the
  one pseudo-case and so adds information rather than removing a prior. The
  per-row figures are on the fitted object as `$smoothing_influence`. Nothing is
  reported once covariates or a grouping variable predict the transitions, since
  those are then fitted by a multinomial logit that the prior never enters —
  which `?fit_lta` now says of `smoothing` generally.

* **The transition prior's size and shape are now documented and sourced.** Both
  were choices and neither was written down: the mass is one pseudo-case per
  origin row rather than per cell, which is the prior Chung, Lanza and Loken
  (2008) use for this model, and it is spread evenly over the destinations
  rather than toward their marginal, because a rare origin row shrunk toward the
  destination marginal would assert that everyone in it moves to the prevalent
  status (Fienberg & Holland, 1973). No default changed.

* **Three small fixes.** A single collapsed pair of latent statuses is now
  reported as "Latent class 1 and 2" rather than "Latent classes"; the
  `@references` blocks of `classification_diagnostics()`, `absolute_fit()` and
  `bivariate_residuals()` had been opened inside `@examples`, which swallowed
  the last example call into the reference text; and the sample-size-adjusted
  BIC now cites Sclove (1987) for its effective sample size.

* **`fit_lta()` and `fit_rmlca()` now accept two-level indicators that are not
  coded 0/1.** `fit_mixture()` has always recoded them for you; the longitudinal
  models never reached that code and stopped with an error instead, so a
  perfectly ordinary 1/2 coding had to be shifted by hand. They now recode it
  themselves, and they decide the mapping once per item across all occasions
  rather than column by column. That distinction matters: the same item appears
  once per occasion, and a per-column decision would map a level differently at
  two occasions whenever one of them happened to observe only one of the two
  levels — which, with thresholds held equal across time, would silently compare
  different response spaces. The recode is reported, naming the item and which
  value became 1, so it is clear which response the printed probabilities
  describe. Items with three or more levels are still an error under
  `measurement = "binary"`.

* **Corrected standard errors now use the estimator they name on models fitted
  with `variances_equal = TRUE`.** `se = "corrected"` needs the sampling
  variance of the measurement parameters, and the vector it built for that
  treated each class's variance as free even on models that hold them equal. The
  extra directions were ones the fit was never able to move along, so the
  numerical information matrix came back indefinite and the package silently
  substituted the outer-product estimator — on a large and very ordinary family
  of models, while printing a diagnostic that read like a failure. The vector
  now carries one variance per item when the classes share it. Standard errors
  for covariate effects on those models change slightly; expect movement in the
  third decimal rather than changed conclusions. `se = "robust"` and
  `se = "hessian"` never used this path and are unaffected.

* **`bayes_constants$categorical` now reaches categorical distal outcomes.** The
  constant was applied to categorical *indicators* but never to a categorical
  *distal outcome*, whose M-step is a separate engine. A class in which nobody
  gave a particular response therefore had no finite estimate for it and the
  intercept ran off towards minus infinity, printing as a large negative logit
  rather than as a bounded one. The prior now enters that M-step in the same
  pseudo-observations form the indicator M-step uses, so the constant means the
  same thing on both, and `categorical = 0` still recovers plain maximum
  likelihood. Estimates for categorical distal outcomes will change, most
  visibly on classes with an unobserved response category. Models with
  covariates alongside the outcome are unaffected: there is no non-arbitrary
  place in covariate space to put the pseudo-observations, so the prior is
  confined to the no-covariate case.

## New: reporting a growth model

* **`measurement_summary()` now works on `fit_gmm()` and `fit_lcga()` fits.** It
  printed a header and nothing else, and returned `NULL`: the growth parameters
  existed only inside `print()`'s output, which cannot be indexed, joined or put
  in a table. It now prints the growth-factor means, the growth-factor variances
  and covariances, the residual variances and the fitted trajectory, and returns
  them in the same long data frame the other models return — `block`,
  `parameter`, `item`, `category`, `class`, `estimate` — with `parameter` taking
  the values `"growth_mean"`, `"growth_variance"`, `"growth_covariance"`,
  `"growth_regression"`, `"residual_variance"` and `"fitted"`. A parameter held
  equal across classes is repeated once per class rather than reported once, so
  the table joins to anything else indexed by class.

* **`compare_longitudinal()` accepts `model = "gmm"` and `model = "lcga"`.**
  Choosing the number of trajectory classes meant looping `fit_gmm()` by hand
  and assembling the table. For the growth models the default `k_range` starts
  at one class — the ordinary latent growth curve model is the benchmark the
  class solutions have to beat, and reporting it is asked for by name in the
  usual reporting checklists. `k_range` now defaults to `NULL`, resolving to
  `1:4` for the growth models and to `2:4`, as before, for `"lta"` and
  `"rmlca"`.

* **`blrt()` takes a fitted growth model with `from_fit =`.** The bootstrap
  likelihood-ratio test has always handled growth mixtures correctly, but
  reaching it meant naming the emission and building the time design with an
  internal function. Passing the fit reads the data, the design, the random
  effects and the covariance constraints off it, so the null and alternative
  models differ from the fit in the number of classes and in nothing else.

* **`lr_test()` warns on a growth model with a collapsed variance**, as it
  already did for the other continuous models. A degenerate fit's
  log-likelihood is not on the same scale as an admissible one, so the test is
  uninterpretable in either direction.

## Fixed: growth mixture models flag a collapsed variance before it reaches zero

* **`fit_gmm()` now judges a variance against the data rather than against the
  estimation floor.** The check fired only once the M-step had pinned a residual
  variance at 1e-6 or a growth-factor covariance eigenvalue at its 1e-8 clip,
  which is the last stage of a collapse rather than the diagnostic one. A class
  with no within-class variation left, or a residual variance three orders of
  magnitude below the others, was therefore reported as an ordinary solution —
  and, because a collapsed variance inflates the likelihood, it could carry the
  best BIC of a whole model set. Both are now compared with the observed
  variance of the outcome on the same 1% rule the latent profile models use, and
  the growth-factor side is judged on the random effects' contribution to the
  implied variance of the outcome rather than on the covariance's own entries,
  which are on the growth factors' scale. A slope variance at zero under a
  healthy intercept variance is still not flagged: that is the
  `random_effects = "intercept"` model, not a degeneracy.

* **The warning says what to do, and what not to do.** It now prints the
  offending variance next to the occasion variance it is small relative to,
  states that this fit's BIC cannot be compared with a clean fit's, and says
  that raising `n_init` will not fix it and can make it worse — a collapse is a
  property of the specification, not of the search, so more starts means more
  chances to find the spike. The flag is stored on the fit and repeated by
  `print()`, so it is still visible on a fit reloaded months later.

## Fixed: two estimation bugs that cost log-likelihood on every affected fit

* **EM now converges before it stops.** Emissions that L-BFGS refines
  afterwards — the binary and Gaussian families — used a much looser stopping
  rule than the others, on the reasoning that the refinement would climb the
  rest of the way. It does not. On a validated four-class binary LCA the loose
  rule stopped 1.27 log-likelihood units short and the refinement recovered
  0.14 of that. The rule also stopped EM at an absolute change of roughly two
  units on a log-likelihood in the thousands, which is far too coarse to *rank*
  restarts, so it was quietly degrading the multi-start search as well. Every
  emission now uses the tight rule, and the refinement starts from a converged
  fit instead of substituting for one. Fits will change, generally for the
  better, and will take more iterations.

* **`bayes_constants` now reaches the L-BFGS refinement.** The refinement read
  `variances` but had `latent` and `categorical` hard-coded at `1`, so setting
  either to `0` switched the prior off in the M-step and left it on in the
  polish. The two stages then optimised different objectives and the polish
  pulled the fit off the maximum-likelihood optimum it had been asked for. The
  documented escape hatch for reproducing an unregularized reference analysis
  now works. The analytical gradient is verified against a finite-difference one
  to 1e-10 for complete data and under FIML, at three prior settings.

* **The L-BFGS refinement no longer runs when class membership is modelled by a
  regression.** Its parameterisation packs a single pooled vector of class
  weights and has no slot for `P(class | covariates)`, so with `predictors` or a
  `group` on the prevalences it was maximising a different model from the one
  being fitted and then writing its measurement parameters back. The damage was
  large and had been invisible: a multiple-group configural model is *separable*,
  so its log-likelihood must equal the sum of the per-group fits, and the polish
  left it 1.10 units short. Removing it closes that to 0.0005 and takes the
  number of restarts reaching the best solution from 1 of 6 to 4 of 6 — it had
  been perturbing every restart away from the optimum, not just the winner.

Together these close a 3.5-unit gap against an external reference on a
three-group latent class model — 5.1 units on the configural model — both now
matched to within 0.001.

## New: binary indicators are recoded to 0/1 for you

* **`measurement = "binary"` now accepts any two-valued item.** A two-level
  factor or character, a logical, and any numeric pair — 1/2, 2/5, whatever the
  source data used — are converted on the way in, and the fit is identical to
  the one hand-coded 0/1 data would have given. Only the mapping is announced,
  once, and `measurement_summary()` then names the level each probability
  belongs to, so "0.87" is never ambiguous. Data already in 0/1 is untouched and
  silent.

* **Two silent-wrong-answer bugs closed on the way.** The old `{0, 1}` check was
  reachable only through a single-string `measurement`, so a 1/2-coded item
  inside a mixed `list()` specification reached the Bernoulli likelihood
  unchecked and returned a wrong log-likelihood with no error. And a *categorical*
  item coded from 0 indexed the previous item's last category rather than its
  own — no error, wrong likelihood. Both are now caught, the second with a
  message saying to add 1.

* An item with three or more values is still refused, but the message now points
  at `"categorical"` or `"continuous"` instead of asking for 0/1 coding that
  would not have helped.

## Renamed: `longitudinal_lrt()` is now `lr_test()`

* The test was never specific to longitudinal models. It takes any two nested
  fits and differences their log-likelihoods and parameter counts, and most uses
  of it — including the multiple-group measurement-invariance test — are
  cross-sectional. `longitudinal_lrt()` still works and warns.

## Changed: what the collapsed-variance warning tells you to do

* The warning used to lead with `bayes_constants = list(variances = 5)`. That
  number was calibrated on one dataset and does not transfer, because the
  constant is divided among the classes: the same value is a different amount of
  prior at every `n_classes`. The warning now leads with
  `variances_equal = TRUE`, which bounds the likelihood so the degeneracy cannot
  arise, then fewer classes, and only then the prior — stated as roughly one
  artificial observation per class and printed as the concrete number for the
  model in hand.

* It also now says two things it should have said before: raising `n_init` is
  **not** a remedy and can make matters worse, since the likelihood really is
  unbounded in that direction and a longer search finds a taller spike; and a
  flagged fit's BIC must not be compared with a clean fit's, because it is
  inflated by the spike. Finally, it asks for a substantive check — that the
  flagged variance is no longer far below the others and that its class mean has
  come off the floor or ceiling — rather than just that the warning stopped.

## New: measurement invariance one parameter at a time

* **New `group_invariant_params` argument on `fit_mixture()`**, for continuous
  indicators. `group_invariant_items` holds whole *items* equal across groups,
  which is the natural constraint when an item has one kind of parameter, as a
  categorical indicator does. A continuous indicator has two — a class mean and
  a variance — and the latent-profile invariance literature routinely frees one
  and holds the other. That model could not be written down before: an item was
  either wholly free across groups or wholly shared.

  `group_invariant_params = "covariances"` frees the class means across groups
  while holding the indicator variances invariant. This is the model
  Olivera-Aguilar and Rikoon (2018) call *unconstrained* and note is the default
  most software fits, and it is the one their invariance test compares against —
  so it is the comparison an applied analysis usually wants, and it is smaller
  and better identified than the fully heterogeneous alternative.
  `group_invariant_params = "means"` is the mirror constraint. The two axes are
  alternatives: pass either `group_invariant_items` or
  `group_invariant_params`, not both.

* **New `variances_equal` argument on `fit_mixture()`**, also for continuous
  indicators: hold each item's variance equal across the classes, so the classes
  differ in location only. This is the homoscedastic latent profile model, and
  the parameterisation several commercial programs estimate by default. It
  applies to an ordinary single-group fit as well, and composes with
  `group_invariant_params` to give a variance shared by the classes but free
  across groups.

  The constrained variance is still stored once per class, so profiles, plots,
  class alignment and the degeneracy check are unchanged; the constraint lives
  in the M-step and in the parameter count.

  A model carrying either constraint is fitted by EM alone. The L-BFGS
  refinement packs one parameter per class-item cell and has no way to express
  an equality across cells, so it would step off the constraint surface; these
  models instead run EM to the tighter stopping rule the unpolished emissions
  already use.

## Improved: the search for a group-varying measurement model

* **`group_effects = "both"` and `"measurement"` no longer rely on random
  starting values alone.** A group-varying measurement model holds one set of
  item parameters per group, tied together only through the class labels, which
  are shared. Random starts give each group its own arbitrary labelling, so
  class 1 in one group and class 1 in another begin as unrelated things and EM
  has no way to bring them into correspondence — and more restarts do not help,
  since every restart draws from the same badly aligned prior. The visible
  symptom was a model that could score *below* the more restricted model nested
  inside it, which is impossible at the optimum and left
  `longitudinal_lrt()` reporting a lower bound instead of a test.

  These models now fit each group on its own data first, permute its classes to
  match the pooled solution, and run one extra restart from there. On the
  `janousch` example the configural fit improves by 15 log-likelihood units and
  the measurement-invariance statistic rises from 251 to 282 on 112 degrees of
  freedom; the substantive conclusion is unchanged, but the number no longer
  depends on the draw. The extra fits make these models roughly a third slower.

  Everything about it degrades gracefully: a group with fewer than
  `2 * n_classes` cases is not fitted on its own, and any sub-fit that fails
  falls back to the pooled parameters, so the worst case is the search as it was
  before.

## Fixed: collapsed class variances in continuous-indicator models

* **A class variance that has collapsed towards zero is now detected and
  reported.** A mixture of normals with freely estimated variances has an
  unbounded likelihood: a class whose variance on one item is driven to zero,
  with its mean on a few cases sharing a value, produces a likelihood that grows
  without limit while describing nothing. Such a solution previously won the
  restart competition, was returned as the best fit, and gave no warning. It now
  raises one naming the item, the class and the group, and reporting whether the
  class mean sits at a data boundary. The flagged cells are stored on the fitted
  object as `$degenerate` and repeated by `print()` and `summary()`, so a saved
  fit still shows the problem months later.

  **This changes results.** A model that previously reported a collapsed
  solution as converged will now warn, and its estimates will differ, because
  the regularisation it is fitted under has changed (see below). The
  `janousch` vignette's measurement-invariance section is the worked example and
  its numbers have changed accordingly.

* **The Gaussian M-step now carries a prior on the variances instead of a hard
  floor.** `gaussian_diag` and `gaussian_diag_nan` previously added `1e-6` to
  each variance after maximising. That is a constant in the units of the data —
  invisible on a five-point scale, enormous on one measured in thousands — and,
  being applied after the fact, did nothing to stop the L-BFGS refinement from
  re-optimising straight through it. The regularisation is now part of the
  objective: a truncated inverse-Wishart prior centred on the item's observed
  marginal variance, expressed as pseudo-observations, applied identically in
  the M-step and in the refinement so the two cannot disagree about what is
  being maximised. On healthy fits the difference is negligible.

* **New `bayes_constants` argument** on `fit_mixture()`, `fit_lta()` and the
  wrappers that delegate to them, exposing the strength of each of the
  estimator's priors — `latent`, `categorical`, `poisson`, `variances` — all
  defaulting to `1`, the value the first three were already hard-coded to.
  Raising `variances` is one of the remedies the collapsed-variance warning
  offers; setting a constant to `0` removes that prior and recovers plain
  maximum likelihood for that block, which is an escape hatch for reproducing an
  unregularized reference analysis rather than a recommended setting.

* **`longitudinal_lrt()` refuses to interpret a degenerate fit.** It warns when
  either input carries the flag — the statistic is meaningless in either
  direction, since a degenerate log-likelihood reflects a spike in the
  likelihood rather than fit to the data. Its existing warning for a negative
  statistic now names both things that can cause one: a failed search on the
  full model, or a degenerate restricted model outscoring it.

* **The `janousch` measurement-invariance comparison was testing the wrong
  restriction.** Measurement invariance is the hypothesis that the item
  parameters are equal across groups, which is `group_effects = "prevalence"`
  (measurement pooled, prevalences free) against `"both"`. The vignette and test
  compared `"measurement"` against `"both"`, which frees the item parameters and
  pools the prevalences — the opposite restriction. Both have been rewritten.

  This supersedes the note in the 0.2.0-era commit that introduced a fixed seed
  for that comparison: the seed was pinning a degenerate configural solution
  rather than curing a search failure, and it has been removed.

## New: `plot()` on a model-selection sweep

`compare_mixtures()` and `compare_longitudinal()` now return an object of class
`mixture_comparison`, and `plot()` on it draws the information criteria against
the number of classes — the elbow plot the fit table was already being read as.
The returned object indexes exactly as the plain list it was, so
`result$fit_table`, `result$models` and `result$best_k` are unchanged.

`indices` takes any subset of `"BIC"`, `"AIC"` and `"SABIC"`; the
log-likelihood and the entropy are deliberately not allowed on that axis, being
on a different scale. `entropy = TRUE` puts relative entropy in a second panel
below on a fixed 0-1 axis rather than on a twin axis, which would invite reading
it as a fit criterion. Each line's minimum is marked, and a K whose maximum was
found by a single random start is drawn hollow. Following Masyn (2013), the plot
is for reading the point of diminishing returns rather than obeying the
minimum — the BIC often keeps falling slowly as classes are added.

## The coefficients and standard errors are now recoverable on the log scale

No printed output changes anywhere in this group. The odds ratios, intervals and
p-values in `print()`, `summary()` and `confint()` are what applied researchers
report and they are shown exactly as before. What changes is that the log-scale
quantities behind them can now be got at.

* **`confint()` no longer rounds inside the object it returns.** It rounded the
  odds ratio and both bounds to three decimals *in the data*, not just for
  display, which destroyed precision in a stored result and put a 0.001 floor
  under any comparison of these numbers against another program's — larger than
  the disagreement such a comparison is usually trying to measure. Values are
  now returned at full precision and rounded only by the print method.
* **New `vcov()` method**, so `sqrt(diag(vcov(fit)))` gives the standard errors
  of the class-membership coefficients. It returns the `(K - 1) * D` matrix over
  the free coefficients, names its rows and columns `"Class k:predictor"`, and
  carries the estimator's name in a `method` attribute as `confint()` already
  does. The class the coefficients are taken against is reported in a
  `ref_class` attribute rather than assumed, since the classes are reordered by
  size after estimation.
* **`coef()` gains `exponentiate`.** The default `TRUE` is exactly the existing
  behaviour. `FALSE` returns the multinomial-logit coefficients themselves. This
  is only convenience over `log(coef(fit))`, which already worked, but it makes
  the log scale discoverable from the documentation.

## Minor improvements

* The collapsed-class-variance warning is shorter, so it is no longer cut off by
  R's default `warning.length` of 1000 bytes — which is what most consoles use,
  and which meant the remedies at the end of the message were the part that
  disappeared. It now names the flagged cells, says the fit is not interpretable
  and its BIC not comparable, and lists the three remedies; the reasoning it
  dropped was already in `?fit_mixture`, which it points at. No default and no
  number changes.

# mixtureEM 0.2.0

## New: stepwise analyses on a fitted model

* **`add_covariates()` and `add_outcome()`.** Once you have chosen an
  unconditional model, relating its classes to external variables no longer
  requires re-specifying (and re-estimating) the whole measurement model:

  ```r
  fit <- fit_mixture(items, n_classes = 3)
  summary(add_covariates(fit, data[, c("age", "sex")]))
  summary(add_outcome(fit, data$bmi))
  ```

  Both verbs run only steps 2–3 of the bias-adjusted three-step approach
  (Bakk et al., 2014; Vermunt, 2010) on the stored step-1
  solution, so the classes you inspected are exactly the classes the
  structural model describes — same solution, same class order, and no risk
  of a re-run landing in a different optimum. `fit_mixture(predictors = )`
  and `fit_mixture(outcome = )` remain available for one-call workflows and
  one-step (simultaneous) estimation.

## New: summary functions return their numbers

* `measurement_summary()` now returns (invisibly) a long-format data frame of
  the item parameters — `item`, `category`, `class`, `estimate` — alongside
  the printed tables, so results can be reused programmatically without
  touching the fitted object's internals.
* `summary()` on a model with a structural part returns (invisibly) a list of
  tidy tables: `$coefficients` (class contrasts with SEs, p-values, odds
  ratios and confidence limits), `$omnibus` (per-covariate Wald tests), and
  `$outcome` (predicted probabilities or class means with their tests).
* New `class_sizes()`: one row per class with the model-estimated proportion,
  the expected count, and the modal-assignment count — the three numbers
  applied papers report.

## Output and plotting improvements

* Long variable and category names no longer break table alignment. Label
  columns size themselves to their content, and names past ~28 characters are
  abbreviated in the display with a key printed under the table (the returned
  data frames always carry the full names).
* Profile and trajectory plots no longer clip long class labels on the right:
  margins are computed from the actual label widths, the legend position
  scales with the x-range, and overly long labels are shortened. The bottom
  margin now adapts to long rotated indicator names as well.
* Factor covariates: levels with no observations (e.g. after subsetting a
  pooled dataset) are dropped automatically instead of producing a degenerate
  `OR = 1.00 [1.00, 1.00]` row, and levels observed fewer than 5 times
  trigger a warning that their coefficients will be unstable.

## Breaking changes

* **Standard errors for covariate effects in 2- and 3-step models have
  changed, and every p-value and confidence interval for a class predictor
  moves with them.** They were computed from the Hessian the M-step stores,
  which measures the curvature of the *Q function* rather than of the step-3
  log-likelihood the estimates actually maximise, and which takes no account
  of the fact that proportional assignment gives each case K weighted
  records. Both errors run the same way, so the reported standard errors were
  too small and the intervals too narrow.

  The new default, `se = "corrected"` on `fit_mixture()`, is the first-order
  corrected estimator of Bakk, Oberski and Vermunt (2014): the step-3
  sandwich plus the variance carried over from step 1, which is estimated and
  not known. `se = "robust"` gives the sandwich alone; `se = "hessian"`
  inverts the step-3 observed information. See `?covariate_se`.

  A 250-replication coverage study on the design of Bakk et al.
  (`data-raw/covariate_se_simulation.R`) says how much this mattered. With
  n = 500 and entropy R² = 0.64, the old nominal 95% intervals covered
  between **54% and 88%** of the time, worst for the strongest coefficient;
  the new default covers 92–98%. With entropy R² near 0.90 the old intervals
  covered 87–95% and the new ones 94–98%, and there the correction term is
  negligible — the sandwich alone does the work. **If you have published
  covariate p-values from this package with small or poorly separated
  classes, they were anti-conservative and are worth re-running.**

  `summary()`, `confint()` and `analytical_wald_test()` now print the name of
  the estimator behind the numbers they show. `correction = "BCH"`, one-step
  fits, and a covariate combined with a distal outcome in one structural
  model keep the uncorrected Hessian, and say so.

## New features

* **`summary()` prints an omnibus Wald test for each class predictor.** For a
  covariate with more than two categories, none of the per-dummy rows tests
  the covariate itself; the new block, printed under the coefficients, gives
  one test per covariate on (K−1) × (levels−1) degrees of freedom and
  inherits whichever variance estimator the model was fitted with.
  `prepare_covariates()` now records which dummy columns belong to one
  variable, and `analytical_wald_test()` matches `term_name` against the
  *variable* rather than by substring. The printed block carries the
  Hauck–Donner caveat: a covariate strong enough to nearly separate a class
  drives its own Wald statistic back towards zero.

* **Mixture latent Markov models, and the mover-stayer restriction.**
  `fit_lta()` gains `n_classes`, which puts latent classes above the chain:
  each class gets its own initial distribution and its own transition
  matrices while they share one measurement model. `mover_stayer = TRUE`
  restricts the last class to the identity transition matrix — a group with
  zero probability of change (Vermunt, 2004). It answers a
  question a single chain cannot: is this one process, or a mobile group
  alongside a stable one?

  `transition_matrix()`, `status_prevalences()` and `plot()` take a `class`
  argument; `print()` and `summary()` report the class sizes and each class's
  chain. Standard errors are **not** available for a mixture over chains, and
  the *unrestricted* mixture is demanding of the data — with a single
  indicator per occasion it is routinely not identified, and classes that
  collapse onto the same chain raise a warning.

* **Absolute and local fit diagnostics for categorical models.**
  `absolute_fit()` compares observed response-pattern frequencies with the
  model-implied ones (likelihood-ratio L², Pearson X², Cressie–Read).
  `bivariate_residuals()` is the local counterpart: the Pearson chi-square of
  each *pair* of items against the model-implied two-way table, divided by
  its degrees of freedom — the conditional-independence assumption showing
  its seams (Oberski et al., 2013). `plot()` draws them
  as a heat table anchored at 1. `classification_table()` yields the
  classification error that the bias-adjusted 3-step estimators exist to
  undo, and `classification_diagnostics()` prints it alongside the AvePP
  matrix. Neither diagnostic enumerates the response-pattern table, so both
  stay usable as indicators are added; the tests re-derive each statistic the
  slow way on small tables and require agreement to machine precision.

* `classification_diagnostics()` now uses the case weights, which matters for
  frequency-weighted pattern files where one row stands for many cases.

* **`fit_gmm()` fits growth mixture models**: class-specific growth curves in
  which cases also vary *within* their class through random growth factors
  (Muthén & Shedden, 1999). Because the outcome is continuous the random
  effects integrate out in closed form, so everything the package already
  offers applies unchanged: `compare_mixtures()` and `blrt()`, predictors and
  distal outcomes through the three-step machinery, survey designs, and FIML
  for unobserved waves. The covariance structure is set by `random_effects`,
  `psi`, `residual`, and `residual_equal`, with conventional defaults.
  Solutions at the edge of the parameter space — a residual variance driven
  to zero, or a singular growth-factor covariance — warn rather than being
  reported as ordinary estimates.

* `fit_gmm()` also takes `growth_predictors`: covariates on the growth
  factors themselves, answering who *within* a trajectory class starts higher
  or grows faster — a different question from `predictors`, which asks who is
  in the class at all. Both can be asked in the same model with
  `n_steps = 1`. Growth-factor covariates must be complete: a case missing
  one has no model-implied mean, and `fit_gmm()` says so rather than deleting
  the case.

* **`fit_lcga()` fits latent class growth analysis**: each class follows its
  own fixed-effect polynomial trajectory (Nagin, 2005), with
  `family = "binomial"`, `"gaussian"`, or `"poisson"` setting the link scale.
  `degree` sets the shape of the curve and `time_scores` the spacing of the
  occasions. `flexmix` is a new `Suggests` dependency, used in the tests as
  independent corroboration on simulated designs.

* `plot()` for `fit_lcga()` and `fit_gmm()` draws the observed data behind
  the estimated curves — without it there is no way to tell a class that
  describes a real subgroup from one that has split a single cloud in half.

* `measurement = "count"` fits a Poisson measurement model: each indicator is
  Poisson-distributed within class with its own rate. Counts work everywhere
  the other indicator types do, and `measurement_summary()` reports a rate
  (lambda) per item and class.

* An emission can ask the EM driver for a higher iteration ceiling and a
  two-stage multi-start search, in which a short first pass ranks the
  restarts and only the survivors are run to convergence. Only the growth
  mixture emission uses them; every other model runs exactly the loop it did
  before.

## Bug fixes

* `fit_lta()` reported a transition as sitting on the zero boundary only when
  it fell below `1e-6`; EM stops when the log-likelihood stops moving, not
  when a cell reaches zero, so practically-zero cells went unflagged. The
  threshold is now `1e-4`.

* **EM now converges for count, polytomous and mixed measurement models.**
  The default stopping rule is loose on purpose for binary and continuous
  indicators, whose fits are polished by an L-BFGS refinement; every other
  measurement model was outside that refinement, so where EM stopped was the
  answer — after five to ten iterations, short of the optimum by enough to
  reverse a BIC comparison between class solutions. Affected models now run
  the iterations they always needed (and are accordingly slower), including
  `blrt()` replicates on those models.

* `fit_mixture()` warns when EM stops at `max_iter` without converging,
  instead of reporting it only through `print()`.

* Cases missing on *every* indicator are deleted before estimation and
  reported by `print()` and in `$missing_data$n_empty_rows`. Previously they
  were retained, which left the parameter estimates correct but inflated *n*
  in BIC/SABIC, assigned every empty case modally to the largest class, and
  dragged relative entropy down. Two of the bundled datasets (`janousch`,
  `yrbs2005`) contain such rows, so their reported sample sizes and entropies
  change accordingly.

## Datasets and vignettes

* The vignette lineup now covers the package's main analyses end to end,
  using only the exported accessor functions: LCA with covariates and a
  distal outcome (`ventura_leon`), latent profile analysis with missing data
  and multiple groups (`janousch`), class enumeration, repeated-measures LCA,
  latent transition analysis with a mover-stayer class, growth mixture
  modelling, and complex-survey LCA (`yrbs2005`).
* The `nlsy79_drinking` dataset and its vignette were removed; the
  repeated-measures LCA vignette now uses a bundled simulated dataset with a
  documented generating script.

# mixtureEM 0.1.0

* Initial release: latent class and profile analysis via mixture modelling,
  including binary, categorical, and continuous measurement models,
  structural models (predictors, distal outcomes), multi-step estimation
  with BCH/ML bias corrections, complex survey designs, missing data
  handling, and longitudinal models (`fit_lta()`, `fit_rmlca()`).
