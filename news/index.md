# Changelog

## mixtureEM (development version)

### Fixed: `absolute_fit()` and `bivariate_residuals()` on a `group_effects = "measurement"` fit

Both statistics used to be computed over the padded J\*Q item matrix a
multiple-group measurement fit stores internally, treating every group’s
copy of each item as a distinct item. That made the implied contingency
table `prod(levels)^Q` cells instead of the `Q * prod(levels)` a
`P(y | group)` model actually implies – silently, since the wrong table
is still a real table and the wrong `df` still produces a p-value, just
the wrong one (an inflated `df` from the padded table understating the
discrepancy between model and data, not a refusal). Both functions now
condition on the group:
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)’s
degrees of freedom are `Q(W - 1) - npar` (which is today’s formula at
`Q = 1`, so an ungrouped fit is unaffected), and
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
returns one item-by-item matrix combining every group’s own pairwise
chi-square rather than a mostly-`NA` matrix duplicated `Q` times. Both
attach the per-group breakdown (`$by_group` / `attr(x, "by_group")`)
alongside the combined statistic. `group_effects = "both"` fits always
carry a covariate structural model for the prevalence effect and were,
and remain, refused outright by the existing conditional-model check –
refit with `group_effects = "measurement"` to check the measurement side
on its own.

### Documentation: the group-dummy Wald test and a `"categorical"` block’s parameter count

[`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`group_effects` argument now warns against testing “is this class the
same size in every group” with
[`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
on the group dummies: that test asks whether a class’s log-odds departs
from the reference group’s, a comparison relative to the group average
rather than a literal constancy test, and the two can disagree.
`group_prevalence_equal` plus
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
is the exact test, and the multiple-group vignette now makes the same
point where it demonstrates that route.

[`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`measurement` argument now documents a real counting limitation: a
single `"categorical"` block uses one `max_val` for the whole block, so
an item with fewer levels than the block’s maximum is still charged for
that maximum, inflating `n_params`, AIC and BIC. The mixed `list`
spelling, splitting items by their own level count, is the remedy.

### Fixed: `n_cores` was undocumented on four entry points

`n_cores` was added to
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md),
[`bootstrap_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md),
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`fit_mixture_internal()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture_internal.md)
when parallel restarts landed, but the `@param` block was only written
for two of the six functions that gained the argument. `R CMD check`
flagged the other four as an `Rd \usage sections` warning. All four now
document what `n_cores` spreads across for that function (bootstrap
replicates or random starts). No code changed and no fitted number
moves.

### Added: partial measurement invariance across groups now works with mixed-type indicators

`group_invariant_items` (and the corresponding partial invariance in
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)/[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md))
used to refuse a mixed (list) `measurement` outright rather than risk
it: the item-to-column bookkeeping and the parameter count both assumed
every item cost the same, which is only true when every item is the same
type. Both now resolve an item to its own sub-model first, so a
partially-invariant model can mix binary and categorical (or continuous)
indicators freely, and a fit’s `n_params` counts each item’s own
parameter cost rather than an average across types. No
previously-reachable fitted number moves.

### Changed: `fit_lta()` now warns when `bayes_constants$latent` is set

[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
has always accepted all four `bayes_constants` names, but the
initial-status and transition priors are governed by the separate
`smoothing` argument, and `latent` was silently ignored. A user
translating another program’s single four-way prior specification could
set it and lose a quarter of it without a word.
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
now warns when `latent` is present in the list, naming `smoothing` as
the argument that actually does the job. No fitted value moves.

### Added: a formula `predictors`, with `data`

[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
gains a `data` argument, and `predictors` can now be a one-sided formula
(`~ age + sex`, or `~ grade * factor(year)` for an interaction)
evaluated against it, the same spelling
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
and
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
already accept. A factor’s dummies and an interaction’s several columns
are recognised as one term by
[`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)’s
omnibus test – exactly the bookkeeping a hand-built
`model.matrix(...)[, -1]` cannot carry, which is what the documentation
used to recommend for a covariate whose effect on class membership is
moderated by another variable.
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
and
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)’s
existing formula support is upgraded the same way: an interaction named
in the formula is now actually expanded and grouped, rather than only
its main-effect variables being extracted.

### Changed: `print()` now points to `ll_knownclass` for a `group=` fit

A fit with a grouping variable has always carried
`metrics$ll_knownclass` and `metrics$n_params_knownclass` – the
log-likelihood and parameter count on the scale software that treats the
group as an observed-without-error class would report – but nothing in
the printed output said so, and the printed `Log-Likelihood` is on the
other scale. [`print()`](https://rdrr.io/r/base/print.html) now adds one
line naming both when they are present. No printed or fitted number
moves otherwise.

### Faster: the E-step no longer walks the sample a case at a time

The E-step normalised its log-likelihoods with a row-wise
[`apply()`](https://rdrr.io/r/base/apply.html), which calls a closure
once per case on every EM iteration of every restart. Profiled on a
13,840-case fit, that one line and its internals were over half the
total runtime. It now uses the package’s own vectorised `logsumexp()`,
which was already sitting in the same source file and already carried a
comment saying [`apply()`](https://rdrr.io/r/base/apply.html) dominates
the runtime wherever it is used. The same substitution was made in the
L-BFGS objective, in the BCH and ML corrections, and in the distal
softmax.

Fits run roughly one and a half to two and a half times faster depending
on the model, and every fitted value is unchanged to the last bit. The
single behavioural difference is a case with zero likelihood under every
class: its normalising constant used to come back `NaN` and is now
`-Inf`, which is the correct value.

### Added: `n_cores`, to fit independent models at the same time

[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md),
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md),
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
and
[`bootstrap_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md)
gain an `n_cores` argument. It defaults to 1, which is the sequential
behaviour these functions have always had. Above 1, the random starts –
or the bootstrap replicates – are spread over that many processes.
[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md),
[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md),
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
pass it through to the engine underneath.

A bootstrap likelihood ratio test at the defaults is over two thousand
model fits, and it is what this helps most: measured three times faster
on four cores.

Results do not depend on `n_cores`. Every random quantity a loop needs –
the starting values of a restart, the synthetic data of a bootstrap
replicate – is drawn in the calling session before any fitting begins,
and a worker only fits what it is handed. The same seed therefore gives
the same answer sequentially and on any number of processes.

Each worker is a fresh R session that has to load the package and
receive a copy of the data, which costs a second or two and some memory.
On a fit that already takes a few seconds, leave `n_cores` at 1; it will
otherwise be slower, not faster.

### Changed: the bootstrap draws `blrt()` takes from a given seed

[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md) now
fixes one seed per replicate before any fitting starts. The replicates
previously shared a single running stream, in which each draw depended
on how many random numbers the previous replicate’s two fits happened to
consume; that is what made the test impossible to reproduce across
process counts. The null distribution obtained from a given seed
therefore differs from the one this function produced before. It is the
same estimator with a different set of draws and the difference is Monte
Carlo error, but a p-value reported from an earlier version will not
reproduce to the digit.

### Removed: the `janousch` dataset and its vignette

The multi-group latent profile example built on `janousch` had grown
into poor teaching material, and both the dataset and its vignette are
gone from the package. The two things it taught – a freed-variance
latent profile analysis with class enumeration, covariates and distal
outcomes, and a multiple-group latent class analysis with a
measurement-invariance test – now live in the two new vignettes
described below. This is a breaking change for any code that calls
`data(janousch)`.

### Added: `liang_ark_sim`, simulated bystander-intervention data

A new bundled dataset with 300 cases, five continuous latent profile
indicators, eight covariates and three continuous distal outcomes. It is
simulated from the published parameters of an existing three-class
latent profile solution and contains no real participant data;
quantities computed on it describe the simulation, not new empirical
findings.

### Added: two vignettes

[`vignette("liang_ark_lpa")`](https://pdvalencia.github.io/mixtureEM/articles/liang_ark_lpa.md)
walks through a full latent profile analysis on `liang_ark_sim`: class
enumeration, what happens when class-specific variances are freed and
collapse, the three-profile solution, covariates, and distal outcomes.
[`vignette("mglca_yrbs")`](https://pdvalencia.github.io/mixtureEM/articles/mglca_yrbs.md)
fits a multiple-group latent class analysis on the bundled `yrbs2005`
data grouped by grade, testing first whether class prevalence differs by
group under a measurement- invariant model and then testing measurement
invariance itself.

### Added: hold a class’s prevalence equal across groups

[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)‘s
new `group_prevalence_equal` argument holds one or more classes’
prevalence to a single shared value across every group while the
remaining classes stay free within each group. This is a restriction the
existing multinomial-logit route (`group_effects = "prevalence"` or
`"both"`) cannot express: a zero coefficient on the group dummies pins a
class to the *reference group’s* prevalence, not to one shared across
every group. The new model is fitted only when a constraint is actually
requested, so an unconstrained fit is unaffected and numerically
equivalent to the existing route. Standard errors for the frozen
prevalences are not produced in this first pass; compare nested models
with
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).

### Added: `class_sizes()` reports each group’s own class sizes

For a fit whose class prevalence varies by group
(`group_effects = "prevalence"` or `"both"`),
[`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md)
now attaches the per-group breakdown as a `"by_group"` attribute,
reconstructed from whichever structural model produced it rather than by
hand.

### Fixed: the bundled `yrbs2005` dataset shipped the raw file, not the analysis sample

`yrbs2005` shipped all 13,917 raw CDC records. Collins & Lanza (2010)’s
analysis sample – the one the vignettes and the documented row count
were already written around – drops the 75 cases missing `grade` as well
as the cases missing on every one of the twelve items, for 13,840 rows,
and the dataset now ships that sample instead of the raw file.
[`vignette("survey_lca")`](https://pdvalencia.github.io/mixtureEM/articles/survey_lca.md)
described a “cases missing on every item” warning firing at fit time;
that prose is gone along with the warning, because those cases are
dropped upstream now rather than surfacing there.

`nrow(yrbs2005)` changes from 13,917 to 13,840. A fit computed on the
old sample is not comparable to one computed on the corrected one.

## mixtureEM 0.3.0

### Fixed: `variances_equal` was ignored by the legacy `n_components` interface

[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
still accepts the older `X`, `n_components`, `Y` and `structural`
argument names, and that spelling is served by its own early return. The
return listed the arguments it forwards explicitly, and
`variances_equal` was not among them; being a named argument rather than
part of `...`, it was not carried through that way either. A
legacy-spelled call that passed `variances_equal` was therefore ignored
without comment, and one that omitted it never picked up the
homoscedastic default a continuous measurement model has carried since
0.2.0. The same analysis written as
`fit_mixture(X, n_components = 2, measurement = "continuous")` and as
`fit_mixture(X, n_classes = 2, measurement = "continuous")` fitted two
different models – free class variances in the first, equal ones in the
second – with nothing in the printed output to say which one you had.

Both spellings now resolve the default and honour an explicit value, and
the check that rejects `variances_equal = TRUE` for indicators it has no
meaning for answers to both. Nothing written with `n_classes` moves, and
neither does any fit of categorical, count or mixed indicators, or any
growth, longitudinal or multiple-group model: those reach the emission
through paths that were never part of this. A continuous-indicator fit
written with `n_components` and no explicit `variances_equal` does move
– it now holds each item’s variance equal across the classes, which is
what the documentation already said it did.

### Standard errors for a continuous distal outcome under BCH were too small

The Wald tests and standard errors reported for a continuous distal
outcome fitted with `correction = "BCH"` were understated, badly enough
to change which effects looked significant and how they ranked against
each other. Reported standard errors get larger; every point estimate is
unchanged.

Two things were wrong. The covariance of the class intercepts and slopes
was built over every row of the data, rather than over the rows with an
observed outcome that the model is actually fitted on – so an outcome
with missing values credited the estimates with information they were
never given. And the standard errors were model-based,
`sigma^2 (U'WU)^-1`, the sandwich estimator having been dropped on the
ground that the bias-corrected fit carries one weighted record per class
per case and a record-level sandwich would count each case several times
over. That is a real problem, but the answer to it is to cluster the
sandwich on the case, not to discard it: the model-based form prices the
fit as if it were ordinary least squares and ignores the information
that inverting the classification table gives up.

The covariance is now built once, alongside the coefficients, from the
same design, weights and mask – so it cannot drift from the estimates it
belongs to – and it is a sandwich clustered on the case. Two checks:
where the class assignment is hard, so that each case contributes a
single record, the estimator reduces exactly to the
heteroskedasticity-consistent (HC0) covariance of the equivalent
one-row-per-case regression; and on a published class-moderation
analysis it now agrees with another program’s standard errors to within
about 2%, against the 2- to 4-fold understatement before.

### `slopes` can now name a subset of covariates for a continuous distal outcome

Previously `slopes` was all-or-nothing: `"pooled"` gave every covariate
one slope shared across classes, `"class_specific"` gave every covariate
its own slope per class. Latent-class moderation analyses routinely need
a mix – the class moderates one set of covariates while the rest are
only adjusted for – and that combination could not be expressed at all.
`slopes` now also accepts a character vector of covariate names (or a
one-sided formula naming them, e.g. `~ loc1 + loc2`): those covariates
get a slope per class, and every other covariate stays pooled. This
applies to
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
and
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`outcome_covariates` path; a categorical outcome still only accepts
`"pooled"` or `"class_specific"`, since doubling the categorical
engine’s surface has no demonstrated need yet.
[`summary()`](https://rdrr.io/r/base/summary.html) gains a third printed
block, “Covariates (Class-specific slopes)”, with one row per covariate
per class and a Wald test of slope equality across classes for each.
Verified two ways: algebraically, the joint design matches
[`lm.wfit()`](https://rdrr.io/r/stats/lmfit.html) on the same expanded
weighted dataset to six decimals, and naming every covariate as
class-specific reproduces `slopes = "class_specific"`’s point estimates
to the same tolerance; and against a real class-moderation analysis with
a mix of moderated and pooled covariates, where every coefficient landed
within about 6% of the reference program’s, all fifteen non-reference
coefficients keeping the same sign.

### `measurement_summary()` gains a `scale` argument for binary indicators

`measurement_summary(fit, scale = "probability")` is unchanged – that
remains the default, and its output is byte-identical to before this
argument existed. Two more scales are now available for a binary
indicator’s item-response probabilities: `scale = "logit"` reports the
same table on the log-odds scale, and `scale = "effect"` reports the
effect-coded parameterisation – an item intercept plus one class
deviation per class, the deviations summing to zero – that several other
programs print by default. This is what makes it possible to place a
mixtureEM measurement model next to such a program’s output at all,
since the two otherwise report different quantities for the same fit;
verified by hand against a reference item’s printed values, matching to
three decimals. A polytomous (more-than-two-category) item is refused
under `scale = "effect"` with a clear error rather than guessed at,
since whether such an item should be coded as ordinal or nominal is a
modelling decision the package does not make on your behalf;
`scale = "logit"` carries no such restriction. The `overall` column,
holding the observed sample marginal, is dropped on both alternative
scales, since a raw proportion has no meaningful transform to either
one.

### Documentation: four clarifications, no behaviour change

[`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
now says outright that the EM convergence rule is fixed and not
user-adjustable, and why: a looser rule was tried and measured to cost
real log-likelihood, and worse, to degrade the multi-start search itself
by ranking candidate starts on numbers too coarse to tell a good
solution from a mediocre one.
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
and
[`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
are cross-referenced – the former has its own, genuinely adjustable
`tol` because a chain mixture needs one, the latter rides
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
fixed rule because it estimates through it.

[`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`group` argument now spells out how to let a covariate’s effect on class
membership vary by group without using `group` at all – build the
interaction directly into `predictors`, e.g.
`model.matrix(~ grade * factor(year))[, -1]` – since class moderation by
a covariate and a `group`-based multiple-group model answer different
questions and are easy to reach for interchangeably by mistake.

[`?bivariate_residuals`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
now says that, for categorical indicators, its statistic agrees closely
with what another program reports under the same name, now that the
expected-count fix above removes the one place they disagreed.
[`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)’s
`n_init` documentation gives a measured runtime figure, so a search of
200 or 1000 starts can be budgeted for rather than guessed at.

### `absolute_fit()` now works with missing data, and `mcar_test()` is new

Until now,
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
simply refused to run once any indicator had a missing value, because
the usual goodness-of-fit table – one cell per possible combination of
answers – cannot be built when different respondents answered different
sets of questions. It now handles this the standard way: cases are
grouped by which items they actually answered, the model is compared to
the data within each such group, and the whole comparison is adjusted
for a saturated (maximally flexible) baseline fit to the same grouping,
so that what remains measures the model rather than the missingness.
This is the “missing at random” (MAR) approach, and it is what gets
printed and returned whenever the data have gaps: the same `g2`, `x2`,
`cressie_read` and `df` as before, now under the weaker assumption, plus
a short block showing the same statistics computed jointly with the
stronger “missing completely at random” (MCAR) assumption for
comparison. A new `dissimilarity` element – roughly, the share of cases
that would need to move to a different response pattern for the model to
fit exactly – is now returned in both the complete-data and the
missing-data case, since a chi-square test on a sparse table so rarely
rejects that a plain descriptive number is often more useful than the
test itself.

`mcar_test(fit)` is new and answers a narrower, separate question:
whether the pattern of missingness itself looks unrelated to the data,
as opposed to whether the model fits. A significant result there does
not mean the model is wrong, and a non-significant one does not mean the
weaker MAR assumption
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
already relies on is safe – the two tests are deliberately independent,
and both are documented as such. Both functions require a plain
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
model with categorical indicators; a model with covariates or continuous
indicators is refused exactly as before, and results on data with no
missing values do not change.

### Fixed: bivariate residuals with missing data flagged the wrong pairs

[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
compares, for every pair of indicators, how often respondents actually
gave each combination of answers against how often the fitted model
expects them to. When one of the two indicators had missing data, the
expected counts were computed from the class sizes for the whole sample,
while the observed counts came only from the respondents who answered
both items — two different groups of people being held to the same
yardstick. The practical effect was that an indicator with a lot of
missing data could make an otherwise unremarkable pair look like the
worst-fitting one in the whole model, simply because of who was missing,
not because the two items were actually related in a way the model
misses. Both counts are now built from the same people: whoever answered
both items in the pair, with the expected counts computed from their own
posterior class membership rather than the sample average. Results are
unchanged when there is no missing data anywhere; with missing data,
some residuals move by an order of magnitude, and the pairs worth
worrying about can change. The printed and returned object now also
reports a single “Total BVR” figure summarizing local dependence across
the whole model, and the bootstrap calibration in `n_reps` respects each
replicate’s missing-data pattern, exact when the missingness is
completely random and an approximation otherwise.

### `outcome_contrasts()`: which classes differ on a distal outcome

The new `outcome_contrasts(fit, ref = NULL)` reports class-vs-class
differences on an outcome attached by
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md),
each with a standard error, a confidence interval and a p-value. With no
`ref` every pair of classes is reported once; naming a class holds it
fixed as the comparison. `adjust` takes `"holm"` or `"bonferroni"` for
the all-pairs table read as a family, and `level` sets the interval.

The omnibus Wald test the summary already printed answers whether the
classes differ at all, which is the question a reviewer asks first and
almost never the one the paper is about. What gets written up is that
the high-risk class scores half a standard deviation above the normative
one and that the two intermediate classes are indistinguishable, and for
a continuous outcome there was no way to get that out of the fit.

Computing it by hand gets it wrong, which is the reason this belongs in
the package. Two class means from one fit are estimated from the same
posteriors — and under the BCH correction from the same inverted
classification-error matrix — so they are correlated, and the standard
error of their difference is not the root of the sum of their squared
standard errors. The contrasts use the full sandwich covariance the fit
already carries, so the joint Wald statistic over the reference
contrasts reproduces the printed omnibus exactly; that identity is
asserted in the tests, for both outcome types, because an indexing error
in a covariance is otherwise silent.

For a continuous outcome,
[`summary()`](https://rdrr.io/r/base/summary.html) now prints the
pairwise differences under the omnibus for up to five classes, and
returns them as `summary(fit)$outcome$contrasts` in every case.
Categorical outcomes already printed a reference-class odds-ratio table
and are unchanged there;
[`outcome_contrasts()`](https://pdvalencia.github.io/mixtureEM/reference/outcome_contrasts.md)
returns the same contrasts in tidy form and adds the pairs that table
never showed. An outcome fitted with `slopes = "class_specific"` is
refused with a message pointing at
[`bootstrap_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md)
and
[`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md):
those models store no covariance between the per-class parameter blocks,
so an exact contrast cannot be formed from them.

### `measurement_summary()` now shows the sample marginal

Every table gains an `Overall` column between the indicator name and the
classes, and the returned data frame gains an `overall` column beside
`estimate`. It is the observed marginal for that item: the weighted
sample proportion beside a probability, the weighted sample mean beside
a mean or a rate, and for a polytomous item the share of cases in each
category.

A conditional number is not readable on its own. A class endorsing an
item at .62 is unremarkable where the sample sits at .60 and is most of
what defines the class where the sample sits at .12, and the table could
not tell those apart — the reader had to go back to the data. Putting
the marginal first in the row makes each class parameter read as a
departure from it.

The benchmark is the *observed* marginal rather than the model-implied
one, because a model-implied column would agree with the class
parameters by construction and so could not act as a check on them. It
uses the case weights where the fit has any, and reads the indicators as
the fit stored them, which is after any binary recode — so the
proportion is of the same level the probability beside it is of. For a
latent transition model the marginal is the one for the occasion being
printed, or, where the measurement model is held equal across occasions,
the one pooled over all of them.

The column is dropped, with a note saying why, for a fit that does not
store its raw indicators and for a multiple-group measurement model,
whose per-group item parameters cannot be matched to the stacked
indicator columns by name; the data frame then carries `NA` there.
Growth models always carry `NA`, a growth factor mean having no sample
marginal to be compared against. No fitted number changes.

### The standardized profile can now be drawn as lines

`plot(fit, type = "line")` draws the same z-scored conditional means as
`type = "bar"`, but as one connected line per class with a zero
reference line. Like the bar chart it needs an all-continuous
measurement model, and it takes the same `scale` argument.

The two are different readings of the same numbers. Bars group by
indicator, so the eye compares classes one indicator at a time; a line
follows a single class across all of them, which is what shows whether
two classes differ in level or in pattern. That is the reading “profile”
names in the applied literature, and until now the only line plot on
offer was the default `type = "profile"`, whose min-max axis has no
meaningful origin and whose shape moves with the sample’s most extreme
observation. For an all-continuous model, prefer `"line"`.

One helper computes the heights for both renderers, so the two can never
disagree about a number.

### Bivariate residuals can now be calibrated by bootstrap

[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
takes `n_reps`, which replaces the statistic’s chi-square reference with
a parametric bootstrap and attaches a matrix of p-values that the print
method shows beside each residual.

The reference distribution is the problem being solved. In the
simulation of Oberski, van Kollenburg and Vermunt (2013), a bivariate
residual referred to chi-square rejected at a nominal five percent in
zero of two hundred samples in *every one* of eight null conditions; its
empirical mean was between 0.25 and 0.36 where the reference has 1. A
statistic that never rejects under a true model also says little under a
false one, so a low bivariate residual is not evidence of good fit — and
the ranking it supports is the least powerful of the three methods those
authors compared. The documentation now says all of this plainly, along
with two limits it had not stated: power falls as the classes separate,
and under a missing-at-random mechanism a large residual is ambiguous
between local dependence and selection bias.

The default `n_reps = 0` is the previous behaviour at the previous cost,
and returns bit-identical residuals. `100` is the recommended working
value and `500` the publication-grade one.

### A standardized profile bar chart

`plot(fit, type = "bar")` draws the grouped bar chart applied
latent-profile papers publish: indicators along the x-axis, one bar per
class within each group, and a zero line separating above-average from
below-average conditional means. It requires an all-continuous
measurement model; `type` defaults to `"profile"`, so the existing
figure is unchanged.

Bar heights are z-scores rather than the profile plot’s min-max scaling,
which is hostage to a single extreme observation and has no meaningful
origin. A z-score is scale-free, so the figure comes out the same
whether or not the indicators were standardized before fitting —
standardizing becomes a display choice rather than a step in preparing
the data. `scale = "within"` divides by the model-implied within-class
standard deviation instead, giving a Cohen’s-d-like reading against
residual rather than total dispersion.

### Step-3 standard errors: which one to compare against another program

The `se` documentation now records that `"corrected"`, the default, is
the statistically right answer, while `"robust"` is the *comparability*
setting. Another program reports the step-3 sandwich alone, so
reproducing its standard errors requires asking for `"robust"`; under
the default a user checking mixtureEM against it sees wider intervals,
and that difference is a difference in estimator — the corrected form
carries step-1 uncertainty the sandwich omits — rather than a bug in
either program. No estimates or standard errors change.

### `measurement` is now required

[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
no longer default `measurement` to `"binary"`. Omitting it is an error
that lists the valid types, shows the mixed-type syntax, and suggests a
type read off your columns:

    `measurement` must be specified. Valid types: "binary", "categorical",
    "continuous", "count".
    Your 8 indicator columns all take two values, so you probably want
      measurement = "binary"
    For mixed types:
      measurement = list(binary = 1:5, continuous = 6:8)

The suggestion is a hint to confirm, not a choice the package makes. The
storage mode of a column does not determine its measurement model: a 1-5
column is a legitimate `"categorical"`, `"continuous"` or `"count"`
indicator, and the class solution differs across the three. Inferring
the type would settle a modelling question by inspecting storage mode
and would make a script’s meaning depend on the data it is run against;
a constant default is the same guess with the data-dependence removed.

This breaks any call that relied on the default. The fix is to add
`measurement = "binary"`, which reproduces the previous behaviour
exactly. `blrt(from_fit = )` is unaffected, since it reads the
specification off the fitted model.

### Continuous indicators default to equal variances across classes

`fit_mixture(measurement = "continuous")` and
`compare_mixtures(measurement = "continuous")` now fit the homoscedastic
latent profile model, holding each item’s variance equal across the
classes. This changes the estimates, the fit indices and the class
solution of any continuous fit that did not set `variances_equal`.
`variances_equal = FALSE` recovers the previous behaviour exactly.

The reason is that the unrestricted normal-mixture likelihood is
unbounded: send a class mean to any single data point and that class’s
variance to zero and the likelihood diverges, so no maximum likelihood
estimate exists and what the EM reports is a local optimum. Holding the
variances equal bounds the likelihood, and the constrained estimator is
consistent (Day, 1969; Hathaway, 1985). Freeing them also invites
classes that describe non-normality in a single population rather than
distinct subgroups (Bauer & Curran, 2003).

The restriction is substantive and testable, and the expectation is that
you fit both and compare. It is not the safe choice but the well-posed
one: the homoscedastic model fails visibly, by splitting a genuinely
heteroscedastic class in two, while the free model fails silently, as a
boundary solution that gets written up as a finding.

Only the two user-facing entry points resolve this default. The growth,
time-block and group-block paths — LCGA, GMM and RMLCA — are unchanged.

### `add_covariates()` and `add_outcome()` accept a formula and `data`

Both functions gain a `data` argument, so the columns can be named
instead of extracted:

``` r

add_covariates(fit, ~ T1age + T1sex + T1SHexp, data = df)
add_outcome(fit, ~ T3NHR_1, data = df)
```

`predictors` and `covariates` also accept a character vector of column
names.
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)’s
formula must name exactly one column, since one call fits one distal
outcome.

This is a matter of typing, and nothing more. The existing calling style
— `add_covariates(fit, df[, c("T1age", "T1sex")])`,
`add_outcome(fit, df$T3NHR_1)` — is unaffected and unchanged, with no
deprecation and no message: passing a computed vector such as `scale(y)`
is often the right thing to do, because the variable is not a column of
anything. Both forms meet at the same code as soon as the columns are in
hand, and every estimate is identical either way.

### Two defaults change: the class-membership prior now reaches step 3

`bayes_constants$latent` is documented as the Dirichlet prior on the
class probabilities. It applied when those probabilities were estimated
as K-1 free weights, and silently stopped applying the moment covariates
entered and they became a regression instead. That leak is now closed:
the prior applies to the class-membership regression too, written as
fractional pseudo-data — one row per class per unique covariate pattern,
weight `latent / (K * U)`, adding `latent / K` cases to each class.

Two consequences, both deliberate:

- **Covariate coefficients move.** Every one of them shrinks slightly
  towards zero, because the prior makes the class sizes slightly more
  equal.
- **Their standard errors move.** They shrink, because unlike the
  “ghost” observation the prior replaces, these rows are part of the
  objective being maximised and so enter the information matrix.

`bayes_constants = list(latent = 0)` restores the previous behaviour
exactly, in both the coefficients and the standard errors; the ghost
observation that guards against complete separation is kept in that
case, as before.

### The three-step corrections take an `assignment` argument

[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md),
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md)
and
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
gain `assignment = c("proportional", "modal")`: how step 1’s posteriors
are turned into the assigned-class variable whose classification error
the BCH and ML corrections invert. **The default does not change** —
`"proportional"` follows Bakk, Tekle and Vermunt (2013), who found it at
least as accurate as modal assignment across 54 simulation conditions
and clearly better when the classes are poorly separated. Use
`assignment = "modal"` when reproducing an analysis whose classes were
assigned that way. The rule in force is stored on the fit and printed
next to the correction, so a saved model still says which one produced
it.

[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
also now documents two things that matter when comparing coefficients
with a published set: a case missing a predictor is retained and
completed under the class-invariant Gaussian marginal rather than
listwise deleted, so the analysed N can differ; and `se = "corrected"`
carries the step-1 uncertainty where `se = "robust"` reports only the
step-3 sampling variability.

### New: `class_assignments()`

The per-case classification now has an accessor, so reaching it no
longer means writing `max.col(fit$log_resp)` by hand. `type = "modal"`
gives the assigned class, `"posterior"` the full matrix, and `"both"` a
data frame carrying the assignment together with its probability — a
per-case classification certainty. It works on `mixture_model`, the
growth models, and `lta_model`, where the assignment is of latent status
and an `occasion` argument picks one out.

Its documentation carries the warning that is the reason it exists: do
not use the returned class as though it were an observed variable in a
subsequent regression, ANOVA or t-test. Use
[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
and
[`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md),
which correct for the classification error that discards.

### `print()` now shows the full set of fit indices

[`print()`](https://rdrr.io/r/base/print.html) on a fitted model showed
the log-likelihood and relative entropy; reading its BIC meant reaching
into `fit$metrics` or running a one-model
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md).
It now prints `Log-Likelihood`, `Parameters`, `AIC`, `BIC`, `SABIC` and
`Rel. Entropy` — exactly the columns
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
tabulates, so printing one model and comparing a range of K can never
show two different sets of numbers for the same fit.
[`print()`](https://rdrr.io/r/base/print.html) on a latent transition
model gains the two indices it was missing.

On a three-step fit the criteria are read off the same set of metrics as
the log-likelihood above them, never a mixture of the two. On a fit
whose variances collapsed, the BIC line says so where it appears, since
that number is inflated by the spike and is not comparable with a clean
fit’s.

### The `n_init` advice now scales with the search that actually ran

The “refit with `n_init = 100`” advice fired whenever the maximum was
found by a single start, whatever `n_init` had been — so a user who ran
`n_init = 200` was told to refit with 100. The advice now reads both of
the counts the fit carries, the restarts requested and the restarts run
out to convergence, and says one of three things.

Below 100 requested it is unchanged. Above 100 requested but with fewer
than 100 run out to convergence — which is the ordinary case on a staged
search, where only the most promising restarts are refined — it says
that the maximum failed to replicate among the restarts that were run
out, that this is a thinner test than the requested count makes it
sound, and that `n_init` should be raised further before anything is
read into it. Only when 100 or more restarts reached convergence does it
raise the specification: how well separated the classes are, how heavily
parameterised the within-class structure is, and whether there are more
classes than the data support.

The strongest reading is now hedged where it was asserted. The number of
random starts a mixture needs grows with the number of classes, the
number of free parameters and how poorly the classes separate, so the
message no longer claims that more starts are unlikely to help, and no
longer ranks over-extraction ahead of the other causes. This applies to
[`print()`](https://rdrr.io/r/base/print.html), the warning,
[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md),
which share one helper. No number changes.

### Diagnostics now say what to change, and by how much

No default changed anywhere in this group, so every existing fit returns
the same numbers. What changed is what the package tells you about them.

- **An unreplicated maximum is now a warning, not a line in
  [`print()`](https://rdrr.io/r/base/print.html).** The most informative
  local-maximum signal there is — “the best solution was found by 1 of
  20 starts” — was a [`cat()`](https://rdrr.io/r/base/cat.html) line,
  invisible to anyone working from
  [`summary()`](https://rdrr.io/r/base/summary.html) or from the
  coefficients. It is now a warning that names the remedy: refit with
  `n_init = 100`, and if the maximum still does not replicate, read that
  as a problem with the specification rather than with the search. It
  fires only from ten requested starts upward, below which a lone
  replication carries no information, and it stays silent on a fit that
  has already been flagged for a collapsed variance or a growth-factor
  boundary — those warnings say that raising `n_init` can make matters
  worse, and two warnings must not give opposite advice about the same
  argument.

- **The number of restarts is now reported honestly on the staged
  searches.** A growth mixture model at `psi = "equal"`, or a latent
  transition model with more than one class, ranks its restarts on a
  short pass and carries only three of them to convergence. The report
  counted the survivors, so a 50-start search announced itself as a
  10-start one;
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  kept no counts at all. Both numbers are now carried, and printed as
  “found by 1 of 3 starts that ran to convergence (of 50 requested)”.

- **[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
  and
  [`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
  gain an `Unreplicated` column**, with a line after the table naming
  the class counts to refit before reporting. The per-model warning is
  suppressed inside
  [`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md),
  which would otherwise raise it once per K before the table it belongs
  next to had been printed.

- **The non-convergence warning names a new `max_iter`** — double the
  one that failed — instead of saying “a larger `max_iter`”, and
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  now issues it at all. It has its own EM driver, so it reported
  non-convergence only through
  [`print()`](https://rdrr.io/r/base/print.html).

- **[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
  warns when 100 draws cannot resolve the decision.** A bootstrap
  p-value can only take the values 1/(B+1), 2/(B+1), …, so at 100 draws
  it cannot separate .04 from .06; when the result lands within one step
  of .05 the test now says so and names `n_reps = 999`. It also counts
  draws where the larger model fitted worse than the model nested inside
  it — a symptom of a replicate search that stopped short — reports the
  count as `n_negative`, and recommends `n_init_boot = 50`.

- **The boundary-probability note now offers a way forward**, rather
  than stating the problem and stopping: read it substantively, since an
  item every member of a class answers identically is often the finding,
  or refit with a stronger `bayes_constants = list(categorical = ...)`
  if that parameter needs an interpretable standard error.

- **Every default the estimator chose is now traceable to a source.**
  [`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md)
  gains three sections and a rewritten one: why `n_init = 20` is a floor
  and what the published replication rates actually say about a count of
  1 (including the correction that a 3–10% band puts 1 of 20 *inside*
  it, which is an argument for more starts rather than a verdict); how
  long EM runs and why the doubling escalation; the bootstrap test, its
  p-value formula, and what 100 draws can and cannot resolve; and a
  table of every warning the package raises with the action for each.
  Where a number cannot be sourced — the `1e-4` EM tolerance,
  `n_init_boot = 10` — the vignette says so rather than dressing it up.

- **The growth-model help files now carry the applied reporting
  conventions, with their sources, and no behaviour changed.**
  [`?fit_gmm`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  names the three specification levels applied papers use, says which of
  them `psi` and `residual_equal` correspond to and which this package
  cannot fit, and states the cost of the `psi = "equal"` default — the
  same constraint that buys stability can buy an extra class that is an
  artefact of it.
  [`?fit_lcga`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  says why an LCGA’s information criteria can improve monotonically with
  K, and why that is a symptom rather than a result. Both now say how
  many occasions each polynomial degree needs.
  [`?class_sizes`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md)
  gives the two published small-class conventions and states plainly
  that the package enforces neither, and the two comparison functions
  document how to read the `Entropy` column — anchors, and the fact that
  entropy is not evidence for the number of classes.

- **The `categorical` and `latent` priors now carry their evidence
  too.** The documentation justified the `variances` prior and left the
  other two as bare defaults. Both the strength of one added observation
  and the decision to spread it in agreement with each item’s observed
  marginal — rather than uniformly over the cells — come from Galindo
  Garre and Vermunt (2006), whose simulation finds that form the most
  accurate of those studied and shows why a uniform spread degrades as
  the number of items grows. No value changed; the package already
  implemented the prior their results favour.

- **[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)’s
  `bayes_constants` now reaches the measurement model.** The two priors
  this model uses are documented as a division of labour — `smoothing`
  for the status prevalences and the transition matrices,
  `bayes_constants` for the measurement model — and the code did not
  implement it. `smoothing` was passed into the measurement model’s
  M-step as well, where it took precedence, so
  `bayes_constants = list(categorical = ...)` had no effect at all and
  `smoothing = 0` quietly removed the measurement prior too, returning
  the item-response probabilities of exactly 0 and 1 that prior exists
  to prevent. Fits using the defaults are unchanged to every digit, both
  arguments being 1.

- **[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  now reports how much of each transition row comes from the prior
  rather than from the data**, and says so when it exceeds five
  percentage points on any row, naming the row and the size of the
  effect. An origin status that few cases occupy is a row the one
  pseudo-case of smoothing carries a visible share of, and the share has
  an exact form rather than needing to be estimated. It is a reading
  caution, not a verdict on the fit: those transitions should be
  reported as indicative. The remedy it points at first is
  `transition_invariance = "full"`, which puts every occasion’s cases
  behind the one pseudo-case and so adds information rather than
  removing a prior. The per-row figures are on the fitted object as
  `$smoothing_influence`. Nothing is reported once covariates or a
  grouping variable predict the transitions, since those are then fitted
  by a multinomial logit that the prior never enters — which
  [`?fit_lta`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  now says of `smoothing` generally.

- **The transition prior’s size and shape are now documented and
  sourced.** Both were choices and neither was written down: the mass is
  one pseudo-case per origin row rather than per cell, which is the
  prior Chung, Lanza and Loken

  2008. use for this model, and it is spread evenly over the
        destinations rather than toward their marginal, because a rare
        origin row shrunk toward the destination marginal would assert
        that everyone in it moves to the prevalent status (Fienberg &
        Holland, 1973). No default changed.

- **Three small fixes.** A single collapsed pair of latent statuses is
  now reported as “Latent class 1 and 2” rather than “Latent classes”;
  the `@references` blocks of
  [`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md),
  [`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
  and
  [`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
  had been opened inside `@examples`, which swallowed the last example
  call into the reference text; and the sample-size-adjusted BIC now
  cites Sclove (1987) for its effective sample size.

- **[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  and
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)
  now accept two-level indicators that are not coded 0/1.**
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  has always recoded them for you; the longitudinal models never reached
  that code and stopped with an error instead, so a perfectly ordinary
  1/2 coding had to be shifted by hand. They now recode it themselves,
  and they decide the mapping once per item across all occasions rather
  than column by column. That distinction matters: the same item appears
  once per occasion, and a per-column decision would map a level
  differently at two occasions whenever one of them happened to observe
  only one of the two levels — which, with thresholds held equal across
  time, would silently compare different response spaces. The recode is
  reported, naming the item and which value became 1, so it is clear
  which response the printed probabilities describe. Items with three or
  more levels are still an error under `measurement = "binary"`.

- **Corrected standard errors now use the estimator they name on models
  fitted with `variances_equal = TRUE`.** `se = "corrected"` needs the
  sampling variance of the measurement parameters, and the vector it
  built for that treated each class’s variance as free even on models
  that hold them equal. The extra directions were ones the fit was never
  able to move along, so the numerical information matrix came back
  indefinite and the package silently substituted the outer-product
  estimator — on a large and very ordinary family of models, while
  printing a diagnostic that read like a failure. The vector now carries
  one variance per item when the classes share it. Standard errors for
  covariate effects on those models change slightly; expect movement in
  the third decimal rather than changed conclusions. `se = "robust"` and
  `se = "hessian"` never used this path and are unaffected.

- **`bayes_constants$categorical` now reaches categorical distal
  outcomes.** The constant was applied to categorical *indicators* but
  never to a categorical *distal outcome*, whose M-step is a separate
  engine. A class in which nobody gave a particular response therefore
  had no finite estimate for it and the intercept ran off towards minus
  infinity, printing as a large negative logit rather than as a bounded
  one. The prior now enters that M-step in the same pseudo-observations
  form the indicator M-step uses, so the constant means the same thing
  on both, and `categorical = 0` still recovers plain maximum
  likelihood. Estimates for categorical distal outcomes will change,
  most visibly on classes with an unobserved response category. Models
  with covariates alongside the outcome are unaffected: there is no
  non-arbitrary place in covariate space to put the pseudo-observations,
  so the prior is confined to the no-covariate case.

### New: reporting a growth model

- **[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  now works on
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  and
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  fits.** It printed a header and nothing else, and returned `NULL`: the
  growth parameters existed only inside
  [`print()`](https://rdrr.io/r/base/print.html)’s output, which cannot
  be indexed, joined or put in a table. It now prints the growth-factor
  means, the growth-factor variances and covariances, the residual
  variances and the fitted trajectory, and returns them in the same long
  data frame the other models return — `block`, `parameter`, `item`,
  `category`, `class`, `estimate` — with `parameter` taking the values
  `"growth_mean"`, `"growth_variance"`, `"growth_covariance"`,
  `"growth_regression"`, `"residual_variance"` and `"fitted"`. A
  parameter held equal across classes is repeated once per class rather
  than reported once, so the table joins to anything else indexed by
  class.

- **[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
  accepts `model = "gmm"` and `model = "lcga"`.** Choosing the number of
  trajectory classes meant looping
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  by hand and assembling the table. For the growth models the default
  `k_range` starts at one class — the ordinary latent growth curve model
  is the benchmark the class solutions have to beat, and reporting it is
  asked for by name in the usual reporting checklists. `k_range` now
  defaults to `NULL`, resolving to `1:4` for the growth models and to
  `2:4`, as before, for `"lta"` and `"rmlca"`.

- **[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
  takes a fitted growth model with `from_fit =`.** The bootstrap
  likelihood-ratio test has always handled growth mixtures correctly,
  but reaching it meant naming the emission and building the time design
  with an internal function. Passing the fit reads the data, the design,
  the random effects and the covariance constraints off it, so the null
  and alternative models differ from the fit in the number of classes
  and in nothing else.

- **[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  warns on a growth model with a collapsed variance**, as it already did
  for the other continuous models. A degenerate fit’s log-likelihood is
  not on the same scale as an admissible one, so the test is
  uninterpretable in either direction.

### Fixed: growth mixture models flag a collapsed variance before it reaches zero

- **[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  now judges a variance against the data rather than against the
  estimation floor.** The check fired only once the M-step had pinned a
  residual variance at 1e-6 or a growth-factor covariance eigenvalue at
  its 1e-8 clip, which is the last stage of a collapse rather than the
  diagnostic one. A class with no within-class variation left, or a
  residual variance three orders of magnitude below the others, was
  therefore reported as an ordinary solution — and, because a collapsed
  variance inflates the likelihood, it could carry the best BIC of a
  whole model set. Both are now compared with the observed variance of
  the outcome on the same 1% rule the latent profile models use, and the
  growth-factor side is judged on the random effects’ contribution to
  the implied variance of the outcome rather than on the covariance’s
  own entries, which are on the growth factors’ scale. A slope variance
  at zero under a healthy intercept variance is still not flagged: that
  is the `random_effects = "intercept"` model, not a degeneracy.

- **The warning says what to do, and what not to do.** It now prints the
  offending variance next to the occasion variance it is small relative
  to, states that this fit’s BIC cannot be compared with a clean fit’s,
  and says that raising `n_init` will not fix it and can make it worse —
  a collapse is a property of the specification, not of the search, so
  more starts means more chances to find the spike. The flag is stored
  on the fit and repeated by
  [`print()`](https://rdrr.io/r/base/print.html), so it is still visible
  on a fit reloaded months later.

### Fixed: two estimation bugs that cost log-likelihood on every affected fit

- **EM now converges before it stops.** Emissions that L-BFGS refines
  afterwards — the binary and Gaussian families — used a much looser
  stopping rule than the others, on the reasoning that the refinement
  would climb the rest of the way. It does not. On a validated
  four-class binary LCA the loose rule stopped 1.27 log-likelihood units
  short and the refinement recovered 0.14 of that. The rule also stopped
  EM at an absolute change of roughly two units on a log-likelihood in
  the thousands, which is far too coarse to *rank* restarts, so it was
  quietly degrading the multi-start search as well. Every emission now
  uses the tight rule, and the refinement starts from a converged fit
  instead of substituting for one. Fits will change, generally for the
  better, and will take more iterations.

- **`bayes_constants` now reaches the L-BFGS refinement.** The
  refinement read `variances` but had `latent` and `categorical`
  hard-coded at `1`, so setting either to `0` switched the prior off in
  the M-step and left it on in the polish. The two stages then optimised
  different objectives and the polish pulled the fit off the
  maximum-likelihood optimum it had been asked for. The documented
  escape hatch for reproducing an unregularized reference analysis now
  works. The analytical gradient is verified against a finite-difference
  one to 1e-10 for complete data and under FIML, at three prior
  settings.

- **The L-BFGS refinement no longer runs when class membership is
  modelled by a regression.** Its parameterisation packs a single pooled
  vector of class weights and has no slot for `P(class | covariates)`,
  so with `predictors` or a `group` on the prevalences it was maximising
  a different model from the one being fitted and then writing its
  measurement parameters back. The damage was large and had been
  invisible: a multiple-group configural model is *separable*, so its
  log-likelihood must equal the sum of the per-group fits, and the
  polish left it 1.10 units short. Removing it closes that to 0.0005 and
  takes the number of restarts reaching the best solution from 1 of 6 to
  4 of 6 — it had been perturbing every restart away from the optimum,
  not just the winner.

Together these close a 3.5-unit gap against an external reference on a
three-group latent class model — 5.1 units on the configural model —
both now matched to within 0.001.

### New: binary indicators are recoded to 0/1 for you

- **`measurement = "binary"` now accepts any two-valued item.** A
  two-level factor or character, a logical, and any numeric pair — 1/2,
  2/5, whatever the source data used — are converted on the way in, and
  the fit is identical to the one hand-coded 0/1 data would have given.
  Only the mapping is announced, once, and
  [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  then names the level each probability belongs to, so “0.87” is never
  ambiguous. Data already in 0/1 is untouched and silent.

- **Two silent-wrong-answer bugs closed on the way.** The old `{0, 1}`
  check was reachable only through a single-string `measurement`, so a
  1/2-coded item inside a mixed
  [`list()`](https://rdrr.io/r/base/list.html) specification reached the
  Bernoulli likelihood unchecked and returned a wrong log-likelihood
  with no error. And a *categorical* item coded from 0 indexed the
  previous item’s last category rather than its own — no error, wrong
  likelihood. Both are now caught, the second with a message saying to
  add 1.

- An item with three or more values is still refused, but the message
  now points at `"categorical"` or `"continuous"` instead of asking for
  0/1 coding that would not have helped.

### Renamed: `longitudinal_lrt()` is now `lr_test()`

- The test was never specific to longitudinal models. It takes any two
  nested fits and differences their log-likelihoods and parameter
  counts, and most uses of it — including the multiple-group
  measurement-invariance test — are cross-sectional.
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  still works and warns.

### Changed: what the collapsed-variance warning tells you to do

- The warning used to lead with `bayes_constants = list(variances = 5)`.
  That number was calibrated on one dataset and does not transfer,
  because the constant is divided among the classes: the same value is a
  different amount of prior at every `n_classes`. The warning now leads
  with `variances_equal = TRUE`, which bounds the likelihood so the
  degeneracy cannot arise, then fewer classes, and only then the prior —
  stated as roughly one artificial observation per class and printed as
  the concrete number for the model in hand.

- It also now says two things it should have said before: raising
  `n_init` is **not** a remedy and can make matters worse, since the
  likelihood really is unbounded in that direction and a longer search
  finds a taller spike; and a flagged fit’s BIC must not be compared
  with a clean fit’s, because it is inflated by the spike. Finally, it
  asks for a substantive check — that the flagged variance is no longer
  far below the others and that its class mean has come off the floor or
  ceiling — rather than just that the warning stopped.

### New: measurement invariance one parameter at a time

- **New `group_invariant_params` argument on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)**,
  for continuous indicators. `group_invariant_items` holds whole *items*
  equal across groups, which is the natural constraint when an item has
  one kind of parameter, as a categorical indicator does. A continuous
  indicator has two — a class mean and a variance — and the
  latent-profile invariance literature routinely frees one and holds the
  other. That model could not be written down before: an item was either
  wholly free across groups or wholly shared.

  `group_invariant_params = "covariances"` frees the class means across
  groups while holding the indicator variances invariant. This is the
  model Olivera-Aguilar and Rikoon (2018) call *unconstrained* and note
  is the default most software fits, and it is the one their invariance
  test compares against — so it is the comparison an applied analysis
  usually wants, and it is smaller and better identified than the fully
  heterogeneous alternative. `group_invariant_params = "means"` is the
  mirror constraint. The two axes are alternatives: pass either
  `group_invariant_items` or `group_invariant_params`, not both.

- **New `variances_equal` argument on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)**,
  also for continuous indicators: hold each item’s variance equal across
  the classes, so the classes differ in location only. This is the
  homoscedastic latent profile model, and the parameterisation several
  commercial programs estimate by default. It applies to an ordinary
  single-group fit as well, and composes with `group_invariant_params`
  to give a variance shared by the classes but free across groups.

  The constrained variance is still stored once per class, so profiles,
  plots, class alignment and the degeneracy check are unchanged; the
  constraint lives in the M-step and in the parameter count.

  A model carrying either constraint is fitted by EM alone. The L-BFGS
  refinement packs one parameter per class-item cell and has no way to
  express an equality across cells, so it would step off the constraint
  surface; these models instead run EM to the tighter stopping rule the
  unpolished emissions already use.

### Improved: the search for a group-varying measurement model

- **`group_effects = "both"` and `"measurement"` no longer rely on
  random starting values alone.** A group-varying measurement model
  holds one set of item parameters per group, tied together only through
  the class labels, which are shared. Random starts give each group its
  own arbitrary labelling, so class 1 in one group and class 1 in
  another begin as unrelated things and EM has no way to bring them into
  correspondence — and more restarts do not help, since every restart
  draws from the same badly aligned prior. The visible symptom was a
  model that could score *below* the more restricted model nested inside
  it, which is impossible at the optimum and left
  [`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  reporting a lower bound instead of a test.

  These models now fit each group on its own data first, permute its
  classes to match the pooled solution, and run one extra restart from
  there. On the `janousch` example the configural fit improves by 15
  log-likelihood units and the measurement-invariance statistic rises
  from 251 to 282 on 112 degrees of freedom; the substantive conclusion
  is unchanged, but the number no longer depends on the draw. The extra
  fits make these models roughly a third slower.

  Everything about it degrades gracefully: a group with fewer than
  `2 * n_classes` cases is not fitted on its own, and any sub-fit that
  fails falls back to the pooled parameters, so the worst case is the
  search as it was before.

### Fixed: collapsed class variances in continuous-indicator models

- **A class variance that has collapsed towards zero is now detected and
  reported.** A mixture of normals with freely estimated variances has
  an unbounded likelihood: a class whose variance on one item is driven
  to zero, with its mean on a few cases sharing a value, produces a
  likelihood that grows without limit while describing nothing. Such a
  solution previously won the restart competition, was returned as the
  best fit, and gave no warning. It now raises one naming the item, the
  class and the group, and reporting whether the class mean sits at a
  data boundary. The flagged cells are stored on the fitted object as
  `$degenerate` and repeated by
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html), so a saved fit
  still shows the problem months later.

  **This changes results.** A model that previously reported a collapsed
  solution as converged will now warn, and its estimates will differ,
  because the regularisation it is fitted under has changed (see below).
  The `janousch` vignette’s measurement-invariance section is the worked
  example and its numbers have changed accordingly.

- **The Gaussian M-step now carries a prior on the variances instead of
  a hard floor.** `gaussian_diag` and `gaussian_diag_nan` previously
  added `1e-6` to each variance after maximising. That is a constant in
  the units of the data — invisible on a five-point scale, enormous on
  one measured in thousands — and, being applied after the fact, did
  nothing to stop the L-BFGS refinement from re-optimising straight
  through it. The regularisation is now part of the objective: a
  truncated inverse-Wishart prior centred on the item’s observed
  marginal variance, expressed as pseudo-observations, applied
  identically in the M-step and in the refinement so the two cannot
  disagree about what is being maximised. On healthy fits the difference
  is negligible.

- **New `bayes_constants` argument** on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  and the wrappers that delegate to them, exposing the strength of each
  of the estimator’s priors — `latent`, `categorical`, `poisson`,
  `variances` — all defaulting to `1`, the value the first three were
  already hard-coded to. Raising `variances` is one of the remedies the
  collapsed-variance warning offers; setting a constant to `0` removes
  that prior and recovers plain maximum likelihood for that block, which
  is an escape hatch for reproducing an unregularized reference analysis
  rather than a recommended setting.

- **[`longitudinal_lrt()`](https://pdvalencia.github.io/mixtureEM/reference/longitudinal_lrt.md)
  refuses to interpret a degenerate fit.** It warns when either input
  carries the flag — the statistic is meaningless in either direction,
  since a degenerate log-likelihood reflects a spike in the likelihood
  rather than fit to the data. Its existing warning for a negative
  statistic now names both things that can cause one: a failed search on
  the full model, or a degenerate restricted model outscoring it.

- **The `janousch` measurement-invariance comparison was testing the
  wrong restriction.** Measurement invariance is the hypothesis that the
  item parameters are equal across groups, which is
  `group_effects = "prevalence"` (measurement pooled, prevalences free)
  against `"both"`. The vignette and test compared `"measurement"`
  against `"both"`, which frees the item parameters and pools the
  prevalences — the opposite restriction. Both have been rewritten.

  This supersedes the note in the 0.2.0-era commit that introduced a
  fixed seed for that comparison: the seed was pinning a degenerate
  configural solution rather than curing a search failure, and it has
  been removed.

### New: `plot()` on a model-selection sweep

[`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
and
[`compare_longitudinal()`](https://pdvalencia.github.io/mixtureEM/reference/compare_longitudinal.md)
now return an object of class `mixture_comparison`, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on it draws the
information criteria against the number of classes — the elbow plot the
fit table was already being read as. The returned object indexes exactly
as the plain list it was, so `result$fit_table`, `result$models` and
`result$best_k` are unchanged.

`indices` takes any subset of `"BIC"`, `"AIC"` and `"SABIC"`; the
log-likelihood and the entropy are deliberately not allowed on that
axis, being on a different scale. `entropy = TRUE` puts relative entropy
in a second panel below on a fixed 0-1 axis rather than on a twin axis,
which would invite reading it as a fit criterion. Each line’s minimum is
marked, and a K whose maximum was found by a single random start is
drawn hollow. Following Masyn (2013), the plot is for reading the point
of diminishing returns rather than obeying the minimum — the BIC often
keeps falling slowly as classes are added.

### The coefficients and standard errors are now recoverable on the log scale

No printed output changes anywhere in this group. The odds ratios,
intervals and p-values in
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) are what applied
researchers report and they are shown exactly as before. What changes is
that the log-scale quantities behind them can now be got at.

- **[`confint()`](https://rdrr.io/r/stats/confint.html) no longer rounds
  inside the object it returns.** It rounded the odds ratio and both
  bounds to three decimals *in the data*, not just for display, which
  destroyed precision in a stored result and put a 0.001 floor under any
  comparison of these numbers against another program’s — larger than
  the disagreement such a comparison is usually trying to measure.
  Values are now returned at full precision and rounded only by the
  print method.
- **New [`vcov()`](https://rdrr.io/r/stats/vcov.html) method**, so
  `sqrt(diag(vcov(fit)))` gives the standard errors of the
  class-membership coefficients. It returns the `(K - 1) * D` matrix
  over the free coefficients, names its rows and columns
  `"Class k:predictor"`, and carries the estimator’s name in a `method`
  attribute as [`confint()`](https://rdrr.io/r/stats/confint.html)
  already does. The class the coefficients are taken against is reported
  in a `ref_class` attribute rather than assumed, since the classes are
  reordered by size after estimation.
- **[`coef()`](https://rdrr.io/r/stats/coef.html) gains
  `exponentiate`.** The default `TRUE` is exactly the existing
  behaviour. `FALSE` returns the multinomial-logit coefficients
  themselves. This is only convenience over `log(coef(fit))`, which
  already worked, but it makes the log scale discoverable from the
  documentation.

### Minor improvements

- The collapsed-class-variance warning is shorter, so it is no longer
  cut off by R’s default `warning.length` of 1000 bytes — which is what
  most consoles use, and which meant the remedies at the end of the
  message were the part that disappeared. It now names the flagged
  cells, says the fit is not interpretable and its BIC not comparable,
  and lists the three remedies; the reasoning it dropped was already in
  [`?fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  which it points at. No default and no number changes.

## mixtureEM 0.2.0

### New: stepwise analyses on a fitted model

- **[`add_covariates()`](https://pdvalencia.github.io/mixtureEM/reference/add_covariates.md)
  and
  [`add_outcome()`](https://pdvalencia.github.io/mixtureEM/reference/add_outcome.md).**
  Once you have chosen an unconditional model, relating its classes to
  external variables no longer requires re-specifying (and
  re-estimating) the whole measurement model:

  ``` r

  fit <- fit_mixture(items, n_classes = 3)
  summary(add_covariates(fit, data[, c("age", "sex")]))
  summary(add_outcome(fit, data$bmi))
  ```

  Both verbs run only steps 2–3 of the bias-adjusted three-step approach
  (Bakk et al., 2014; Vermunt, 2010) on the stored step-1 solution, so
  the classes you inspected are exactly the classes the structural model
  describes — same solution, same class order, and no risk of a re-run
  landing in a different optimum. `fit_mixture(predictors = )` and
  `fit_mixture(outcome = )` remain available for one-call workflows and
  one-step (simultaneous) estimation.

### New: summary functions return their numbers

- [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  now returns (invisibly) a long-format data frame of the item
  parameters — `item`, `category`, `class`, `estimate` — alongside the
  printed tables, so results can be reused programmatically without
  touching the fitted object’s internals.
- [`summary()`](https://rdrr.io/r/base/summary.html) on a model with a
  structural part returns (invisibly) a list of tidy tables:
  `$coefficients` (class contrasts with SEs, p-values, odds ratios and
  confidence limits), `$omnibus` (per-covariate Wald tests), and
  `$outcome` (predicted probabilities or class means with their tests).
- New
  [`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md):
  one row per class with the model-estimated proportion, the expected
  count, and the modal-assignment count — the three numbers applied
  papers report.

### Output and plotting improvements

- Long variable and category names no longer break table alignment.
  Label columns size themselves to their content, and names past ~28
  characters are abbreviated in the display with a key printed under the
  table (the returned data frames always carry the full names).
- Profile and trajectory plots no longer clip long class labels on the
  right: margins are computed from the actual label widths, the legend
  position scales with the x-range, and overly long labels are
  shortened. The bottom margin now adapts to long rotated indicator
  names as well.
- Factor covariates: levels with no observations (e.g. after subsetting
  a pooled dataset) are dropped automatically instead of producing a
  degenerate `OR = 1.00 [1.00, 1.00]` row, and levels observed fewer
  than 5 times trigger a warning that their coefficients will be
  unstable.

### Breaking changes

- **Standard errors for covariate effects in 2- and 3-step models have
  changed, and every p-value and confidence interval for a class
  predictor moves with them.** They were computed from the Hessian the
  M-step stores, which measures the curvature of the *Q function* rather
  than of the step-3 log-likelihood the estimates actually maximise, and
  which takes no account of the fact that proportional assignment gives
  each case K weighted records. Both errors run the same way, so the
  reported standard errors were too small and the intervals too narrow.

  The new default, `se = "corrected"` on
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  is the first-order corrected estimator of Bakk, Oberski and Vermunt
  (2014): the step-3 sandwich plus the variance carried over from step
  1, which is estimated and not known. `se = "robust"` gives the
  sandwich alone; `se = "hessian"` inverts the step-3 observed
  information. See
  [`?covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).

  A 250-replication coverage study on the design of Bakk et al.
  (`data-raw/covariate_se_simulation.R`) says how much this mattered.
  With n = 500 and entropy R² = 0.64, the old nominal 95% intervals
  covered between **54% and 88%** of the time, worst for the strongest
  coefficient; the new default covers 92–98%. With entropy R² near 0.90
  the old intervals covered 87–95% and the new ones 94–98%, and there
  the correction term is negligible — the sandwich alone does the work.
  **If you have published covariate p-values from this package with
  small or poorly separated classes, they were anti-conservative and are
  worth re-running.**

  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
  now print the name of the estimator behind the numbers they show.
  `correction = "BCH"`, one-step fits, and a covariate combined with a
  distal outcome in one structural model keep the uncorrected Hessian,
  and say so.

### New features

- **[`summary()`](https://rdrr.io/r/base/summary.html) prints an omnibus
  Wald test for each class predictor.** For a covariate with more than
  two categories, none of the per-dummy rows tests the covariate itself;
  the new block, printed under the coefficients, gives one test per
  covariate on (K−1) × (levels−1) degrees of freedom and inherits
  whichever variance estimator the model was fitted with.
  `prepare_covariates()` now records which dummy columns belong to one
  variable, and
  [`analytical_wald_test()`](https://pdvalencia.github.io/mixtureEM/reference/analytical_wald_test.md)
  matches `term_name` against the *variable* rather than by substring.
  The printed block carries the Hauck–Donner caveat: a covariate strong
  enough to nearly separate a class drives its own Wald statistic back
  towards zero.

- **Mixture latent Markov models, and the mover-stayer restriction.**
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  gains `n_classes`, which puts latent classes above the chain: each
  class gets its own initial distribution and its own transition
  matrices while they share one measurement model. `mover_stayer = TRUE`
  restricts the last class to the identity transition matrix — a group
  with zero probability of change (Vermunt, 2004). It answers a question
  a single chain cannot: is this one process, or a mobile group
  alongside a stable one?

  [`transition_matrix()`](https://pdvalencia.github.io/mixtureEM/reference/transition_matrix.md),
  [`status_prevalences()`](https://pdvalencia.github.io/mixtureEM/reference/status_prevalences.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) take a
  `class` argument; [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) report the class
  sizes and each class’s chain. Standard errors are **not** available
  for a mixture over chains, and the *unrestricted* mixture is demanding
  of the data — with a single indicator per occasion it is routinely not
  identified, and classes that collapse onto the same chain raise a
  warning.

- **Absolute and local fit diagnostics for categorical models.**
  [`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md)
  compares observed response-pattern frequencies with the model-implied
  ones (likelihood-ratio L², Pearson X², Cressie–Read).
  [`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md)
  is the local counterpart: the Pearson chi-square of each *pair* of
  items against the model-implied two-way table, divided by its degrees
  of freedom — the conditional-independence assumption showing its seams
  (Oberski et al., 2013).
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them as
  a heat table anchored at 1.
  [`classification_table()`](https://pdvalencia.github.io/mixtureEM/reference/classification_table.md)
  yields the classification error that the bias-adjusted 3-step
  estimators exist to undo, and
  [`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md)
  prints it alongside the AvePP matrix. Neither diagnostic enumerates
  the response-pattern table, so both stay usable as indicators are
  added; the tests re-derive each statistic the slow way on small tables
  and require agreement to machine precision.

- [`classification_diagnostics()`](https://pdvalencia.github.io/mixtureEM/reference/classification_diagnostics.md)
  now uses the case weights, which matters for frequency-weighted
  pattern files where one row stands for many cases.

- **[`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  fits growth mixture models**: class-specific growth curves in which
  cases also vary *within* their class through random growth factors
  (Muthén & Shedden, 1999). Because the outcome is continuous the random
  effects integrate out in closed form, so everything the package
  already offers applies unchanged:
  [`compare_mixtures()`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
  and
  [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md),
  predictors and distal outcomes through the three-step machinery,
  survey designs, and FIML for unobserved waves. The covariance
  structure is set by `random_effects`, `psi`, `residual`, and
  `residual_equal`, with conventional defaults. Solutions at the edge of
  the parameter space — a residual variance driven to zero, or a
  singular growth-factor covariance — warn rather than being reported as
  ordinary estimates.

- [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  also takes `growth_predictors`: covariates on the growth factors
  themselves, answering who *within* a trajectory class starts higher or
  grows faster — a different question from `predictors`, which asks who
  is in the class at all. Both can be asked in the same model with
  `n_steps = 1`. Growth-factor covariates must be complete: a case
  missing one has no model-implied mean, and
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  says so rather than deleting the case.

- **[`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  fits latent class growth analysis**: each class follows its own
  fixed-effect polynomial trajectory (Nagin, 2005), with
  `family = "binomial"`, `"gaussian"`, or `"poisson"` setting the link
  scale. `degree` sets the shape of the curve and `time_scores` the
  spacing of the occasions. `flexmix` is a new `Suggests` dependency,
  used in the tests as independent corroboration on simulated designs.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  and
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  draws the observed data behind the estimated curves — without it there
  is no way to tell a class that describes a real subgroup from one that
  has split a single cloud in half.

- `measurement = "count"` fits a Poisson measurement model: each
  indicator is Poisson-distributed within class with its own rate.
  Counts work everywhere the other indicator types do, and
  [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  reports a rate (lambda) per item and class.

- An emission can ask the EM driver for a higher iteration ceiling and a
  two-stage multi-start search, in which a short first pass ranks the
  restarts and only the survivors are run to convergence. Only the
  growth mixture emission uses them; every other model runs exactly the
  loop it did before.

### Bug fixes

- [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  reported a transition as sitting on the zero boundary only when it
  fell below `1e-6`; EM stops when the log-likelihood stops moving, not
  when a cell reaches zero, so practically-zero cells went unflagged.
  The threshold is now `1e-4`.

- **EM now converges for count, polytomous and mixed measurement
  models.** The default stopping rule is loose on purpose for binary and
  continuous indicators, whose fits are polished by an L-BFGS
  refinement; every other measurement model was outside that refinement,
  so where EM stopped was the answer — after five to ten iterations,
  short of the optimum by enough to reverse a BIC comparison between
  class solutions. Affected models now run the iterations they always
  needed (and are accordingly slower), including
  [`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
  replicates on those models.

- [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  warns when EM stops at `max_iter` without converging, instead of
  reporting it only through
  [`print()`](https://rdrr.io/r/base/print.html).

- Cases missing on *every* indicator are deleted before estimation and
  reported by [`print()`](https://rdrr.io/r/base/print.html) and in
  `$missing_data$n_empty_rows`. Previously they were retained, which
  left the parameter estimates correct but inflated *n* in BIC/SABIC,
  assigned every empty case modally to the largest class, and dragged
  relative entropy down. Two of the bundled datasets (`janousch`,
  `yrbs2005`) contain such rows, so their reported sample sizes and
  entropies change accordingly.

### Datasets and vignettes

- The vignette lineup now covers the package’s main analyses end to end,
  using only the exported accessor functions: LCA with covariates and a
  distal outcome (`ventura_leon`), latent profile analysis with missing
  data and multiple groups (`janousch`), class enumeration,
  repeated-measures LCA, latent transition analysis with a mover-stayer
  class, growth mixture modelling, and complex-survey LCA (`yrbs2005`).
- The `nlsy79_drinking` dataset and its vignette were removed; the
  repeated-measures LCA vignette now uses a bundled simulated dataset
  with a documented generating script.

## mixtureEM 0.1.0

- Initial release: latent class and profile analysis via mixture
  modelling, including binary, categorical, and continuous measurement
  models, structural models (predictors, distal outcomes), multi-step
  estimation with BCH/ML bias corrections, complex survey designs,
  missing data handling, and longitudinal models
  ([`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md)).
