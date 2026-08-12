# What mixtureEM Estimates, and Why

``` r

library(mixtureEM)
```

Most of the time you do not need this vignette. The defaults are the
settings the package intends you to use, and the other vignettes show
the analyses.

This one exists for the moments when a number does not look the way you
expected: a log-likelihood that differs from another program’s, a
warning about a collapsed variance, a probability reported as 0.000.
Each of those is a deliberate choice with a literature behind it, and
each is easier to live with once you know which choice it was.

## Mixture likelihoods have more than one maximum

Fitting a mixture model means climbing a surface with several peaks. The
EM algorithm climbs the one it starts nearest, so a single run is a coin
toss rather than an estimate.
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
therefore runs 20 random starts by default and keeps the best, and it
reports how many of them found that best solution:

``` r

set.seed(1)
n <- 400
p <- matrix(c(.9, .85, .8, .15,
              .1, .20, .15, .90), nrow = 2, byrow = TRUE)
z <- sample(1:2, n, TRUE)
X <- matrix(rbinom(n * 4, 1, p[z, ]), n, 4)
colnames(X) <- paste0("q", 1:4)

fit <- fit_mixture(X, n_classes = 2, measurement = "binary", random_state = 1)
fit
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 18 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -909.47
#>   Rel. Entropy   : 0.8316
#>   Best solution  : found by 20 of 20 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 52.06%
#>   Class 2: 47.94%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```

**Why 20.** In Hipp and Bauer’s (2006) simulations of correctly
specified models, the best solution was reached by 44% to 77% of the
converged starts, and “in all replications the highest log-likelihood
occurred in at least 39% of the random starts” (p. 46). At those rates
20 starts should find the maximum somewhere between eight and fifteen
times, which is more than enough to *see* a stable maximum when there is
one. Twenty also sits inside the range other comparisons of
initialisation strategies use: Biernacki, Celeux and Govaert (2003) run
ten repetitions throughout; Shireman, Steinley and Brusco (2017)
recommend thirty configurations.

**How to read the report.** The applied benchmark is a rate, not a
count. Nylund-Gibson and Choi (2018, p. 454) describe a log-likelihood
“replicated 82 times out of the 500 requested, which is more than
sufficient (e.g., some recommend that 3–10% is needed)”. Read carefully,
that band says less about a 20-start search than it looks like it does:
3–10% of 20 is 0.6 to 2 replications, so a count of **1 is inside the
band** and cannot be told apart from luck *at that sample size*. That is
a statement about the resolution of a 20-start search, not a verdict on
the fit — which is exactly why the action is to raise `n_init` rather
than to conclude anything. At 100 starts the band is 3 to 10, and a
count of 1 falls outside it.

**What to do.** Refit with `n_init = 100`. Hipp and Bauer (p. 49) put
the escalation there directly: “simply generating 10 sets of start
values will be insufficient …; instead, at least 50 to 100 sets of
starting values will be needed. This becomes even more necessary with
models containing more classes since they converge less frequently.” If
the maximum still does not replicate at 100 starts, stop raising it. “If
the optimal solution is found only rarely, this may suggest an error in
the model specification” (p. 48) — most often more classes than the data
support. mixtureEM warns about this case rather than leaving it to the
printed line, and the warning names the same number.

**The honest limits.** Twenty is a floor for small, well-separated
models, not a guarantee. Local optima multiply as classes are added and
as separation falls: in Shireman et al.’s simulations the rate of local
solutions rose from 0.62 at two clusters to 0.86 at six, and on one real
dataset of 11,782 cases, 1,000 random starts produced 1,000 *distinct*
solutions — their footnote concludes that “a minimal number of starting
values in order to arrive at an adequate solution should exceed 1,000
for some data types” (p. 284). Nor is agreement among a handful of
starts evidence of anything. A common rule of thumb holds that the best
log-likelihood repeating twice is enough; Hipp and Bauer (p. 49)
explicitly reject it, and in their own case study the best four-class
solution was found by about 3% of converged starts, roughly 1 in 30. Two
agreeing starts is not the bar.

**When two starts count as the same solution.** The replication count
treats restarts whose log-likelihoods are within 0.01 as having found
the same maximum. That is deliberately looser than a numerical-identity
threshold — Shireman et al. separate solutions at a BIC difference of
1e-16 — because the failure it avoids is real and the one it risks is
not: genuinely different optima in these models sit whole units apart,
while four starts landing at -6483.1620, -6483.1613, -6483.1612 and
-6483.1607 are one optimum reported four times. At a tighter threshold
that fit reports “found by 1 of 6”.

Set `random_state` to make the search reproducible.

## How long EM runs, and when it stops

Each start runs until the log-likelihood stops improving by more than
`1e-4`, or until it has taken `max_iter = 1000` iterations, whichever
comes first.

The iteration budget has a direct precedent: Biernacki et al. (2003)
budget 1,000 iterations in their comparison, and the first
recommendation of their discussion is “for a good solution, do not skimp
on the number of iterations” (p. 574). Their exploration moves on a
multiplicative scale — 60, 120, 240, 480, 960 and on up — which is why
the advice when a fit stops at the cap is to **double** `max_iter`
rather than to multiply it by ten. mixtureEM’s non-convergence warning
names the doubled value.

The tolerance is a different kind of number, and it is worth being plain
about it: `1e-4` is a pragmatic calibration against this package’s own
reference fits, not a value taken from the literature. If anything the
literature argues against relying on a rule of that form at all.
Biernacki et al. (p. 568): “We do not use stopping criteria based on the
relative change of the estimates or loglikelihood because the slow
convergence of the EM makes such criteria hazardous.” Their own short
runs stop instead on progress relative to progress already made,
$`(L^q - L^{q-1}) / (L^q - L^0)`$.

That objection is exactly why the hard cap exists alongside the
tolerance rather than instead of it. A slow-converging EM can satisfy an
absolute tolerance while still climbing; the cap makes that visible,
because a fit that hits it says so.

## The bootstrap test and how many draws it needs

[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
compares *k* against *k* − 1 classes by generating data from the smaller
model, refitting both to each generated sample, and locating the
observed likelihood-ratio statistic in the resulting distribution. The
p-value is

``` math
p = \frac{1 + \#\{t^*_r \ge t\}}{B + 1},
```

which is Davison and Hinkley’s (1997) eq. (4.11). The `+1`s are not a
continuity correction: with them, the p-value is uniform on
$`\{1/(B+1), \ldots, 1\}`$ under the null, and “in this sense the Monte
Carlo test is exact” (p. 141).

**Why B = 100 is the default.** Dziak, Lanza and Tan (2014, p. 3):
“Simulation evidence in McLachlan (1987) suggests that B should be at
least 99 to obtain optimal power; we use B = 100 in this paper.” Davison
and Hinkley agree on the floor — “it is advisable to take R to be at
least 99” (p. 143) — and quantify what it costs. Against the full test,
99 replicates retain a power ratio of 0.83 at α = .05 and 0.95 at 999,
so “the loss of power with R = 99 is not serious for α ≥ .05, and R =
999 should generally be safe”. At α = .01 the ratio at 99 draws falls to
0.60.

**Where 100 stops being adequate.** The attainable p-values are 1/(B+1)
apart, so at 100 draws the smallest is 0.0099 and the test cannot
separate .04 from .06. Where the decision turns on that, “the critical
region of a fixed-level test has been randomly displaced” (p. 155).
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
warns when the result lands within one step of .05 and names
`n_reps = 999`.

**`n_init_boot = 10` is a compute compromise, not a recommendation.**
Each draw refits both models, so the replicate search is where the cost
of the test lives, and 10 restarts per replicate is where this package
puts the trade-off. No source endorses that number. Dziak et al. (p. 4,
fn. 3) say the minimum “is not known”, that too few starts under the
alternative “can lead to invalid results”, and that a local maximum
there can make “the calculated likelihood ratio … occasionally be
nonsensical negative value”; they used 50. Lee, Wickrama and O’Neal
(2023, p. 654) make the same point from the applied side: “If a local
solution is detected in either the k-1 class or the k-class model (or
both), the difference between the log-likelihood values will be biased.
In turn, this leads to biased inference statistics (p-values).” That
negative statistic is the diagnostic, and
[`blrt()`](https://pdvalencia.github.io/mixtureEM/reference/blrt.md)
counts it: any negative draw and it warns, naming `n_init_boot = 50`.

**A non-significant BLRT is not evidence that there is no extra class.**
Its power turns on class separation far more than on sample size. In
Tekle, Gudicha and Vermunt’s (2016) Table 2, distinguishing two classes
from three with six indicators needs 25 cases at high separation and 670
at low; for three against four, 1,830. Dziak et al. (p. 6) report power
of 4% at *n* = 100 with unequal classes and poor measurement, against
99% with strong measurement and equal classes. The same test, the same
number of draws.

## When mixtureEM warns you

| Warning | What it means | What to do |
|----|----|----|
| EM did not converge within `max_iter` | The estimates are wherever the algorithm had reached, which need not be a maximum. | Refit with `max_iter` doubled. If doubling does not help, the model is probably weakly identified at this number of classes. |
| The solution was found by 1 start | The maximum may be the best of a small sample of the surface rather than the best there is. | Refit with `n_init = 100`. If it still does not replicate, suspect the specification — most often too many classes. |
| A variance has collapsed | A class has been fitted to a spike rather than to a subgroup. Its likelihood, and so its BIC, is not comparable with a clean fit’s. | Refit with a stronger variance prior, `bayes_constants = list(variances = n_classes)`, or with fewer classes. Do **not** raise `n_init`: more starts is more chances to find the spike. |
| Probabilities at the boundary | An item every member of a class answers identically. The standard error for that parameter is not interpretable. | Read it substantively — this is often the finding — or refit with a stronger `bayes_constants = list(categorical = ...)` if you need the standard error. |
| Latent classes have converged on the same chain | The model is describing one process with the parameters of several. | Refit with fewer classes, or with a restriction such as `mover_stayer = TRUE` that tells them apart. |
| The BLRT p-value is next to .05 | With `n_reps` draws the p-value moves in steps of 1/(B+1) and cannot resolve the decision. | Refit with `n_reps = 999`. |
| Negative bootstrap draws | A replicate’s *k*-class fit landed below the (*k*−1)-class model nested inside it: the replicate search stopped short. | Refit with `n_init_boot = 50`. Until then the p-value is biased. |

## The estimator is a posterior mode, not maximum likelihood

Maximum likelihood in latent class analysis routinely pushes a
conditional probability to exactly 0 or 1. When that happens the
information matrix loses rank, and standard errors, confidence intervals
and significance tests for those parameters stop meaning anything
(Galindo Garre & Vermunt, 2006, sec. 2.1).

mixtureEM avoids that with weak priors that hold estimates strictly
inside the parameter space — the *Bayes constants* approach, controlled
by four named counts:

``` r

bayes_constants = list(latent = 1, categorical = 1, poisson = 1, variances = 1)
```

Each is a count of artificial observations spread over the classes, so
its influence shrinks as the sample grows and is negligible on any class
with real data behind it. All four default to 1.

Three consequences worth knowing:

- **The reported log-likelihood is the log-likelihood**, evaluated at
  the posterior mode — not the log-posterior. It is the same function a
  maximum-likelihood program reports, evaluated at slightly different
  estimates, so it is bounded above by the ML value. The gap grows with
  the number of parameters ML would have placed on the boundary, and it
  keeps AIC and BIC on a defensible scale.
- **Setting a constant to `0` recovers plain maximum likelihood** for
  that block. That is an escape hatch for reproducing an unregularized
  reference analysis, not a recommended setting.
- **[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  flags probabilities that reached the boundary anyway.** Those tell you
  something real — a class defined partly by an item every case in it
  answers the same way — and warn you not to read their standard errors.

## Continuous indicators: the likelihood is unbounded

For continuous indicators the problem is sharper. Put a class mean on a
few identical cases and shrink its variance towards zero, and the
density at those cases grows without limit: the likelihood is
*unbounded*, and the highest-scoring solution can be one that describes
nothing (Day, 1969; Kiefer & Wolfowitz, 1956).

When a fit lands there, mixtureEM says so and names the item and the
class. The important thing about that warning is what it is **not**: it
is not a convergence failure. Raising `n_init` will not fix it and can
make it worse, because a longer search finds a taller spike. For the
same reason a flagged fit’s BIC cannot be compared with a clean one’s —
it is inflated by the spike.

There are three ways out, and they change the model in different ways:

- `variances_equal = TRUE` holds each item’s variance equal across
  classes. The likelihood is then bounded and the problem cannot arise.
  This is the homoscedastic latent-profile model and the constrained
  approach of Hathaway (1985).
- **Fewer classes.** A class describing a spike rather than a subgroup
  usually means the data do not support this many.
- **A stronger prior**, `bayes_constants = list(variances = ...)`. A
  useful starting point is one artificial observation per class — that
  is `variances = n_classes` — increasing it a little at a time if the
  warning persists. The warning prints that number for your model.

Then check the result rather than just that the warning stopped: the
flagged variance should no longer be far below the others, and its class
mean should have come off the floor or ceiling of the response scale.

Two things about that prior are worth knowing, because they are what
make the advice a rule rather than a number to copy.

It is centred on each item’s *own* observed marginal variance, not on a
fixed quantity. That is what lets one constant mean the same thing on a
five-point scale and on an income variable: rescale an indicator and the
prior rescales with it. A penalty pinned to a constant has no such
property, and the penalised-likelihood literature prefers the
data-scaled form for exactly this reason (Chen, Tan & Zhang, 2008,
sec. 4).

And the constant is divided among the classes, so
`variances = n_classes` is the setting that holds it at one artificial
observation per class at *every* number of classes — the way the
literature applies it, the same amount to every component. A bare number
does the opposite: it is a different amount of prior in a three-class
model than in a six-class one. Write it as `n_classes` and it stays
comparable across the models you are enumerating.

Theory settles the *form* of this prior and not its size. Within the
conditions Chen, Tan & Zhang (2008) impose — which any fixed positive
constant meets — a penalty that diverges as a variance approaches zero,
and grows more slowly than the sample, gives a consistent estimator
whatever its constant. That result is proved for *univariate* normal
mixtures; the multivariate case, which is what this package fits, they
leave open (sec. 5). So the choice of constant is all the more a
finite-sample judgement, which is why the default is deliberately weak
and the advice above comes with a check attached.

### A clean fit is not proof that there is no spike

The absence of the warning does not establish that no collapsed solution
exists. The warning can only report on the solution the search actually
reached, and some degenerate solutions sit in basins so narrow that
random starts never land in them. A fit can replicate its maximum on
every one of hundreds of starts, and still not be the highest point on
the likelihood.

This is worth checking whenever an indicator has a floor or a ceiling —
a bounded scale with a visible pile-up of cases at an endpoint. That
pile-up is a set of near-identical values, which is exactly what a
collapsing class latches onto, but it attracts only starts that begin
very close to it. Two cheap checks: refit with
`bayes_constants = list(variances = n_classes)` and see whether the
solution moves at all, and look at whether any class mean is sitting on
a scale endpoint with a variance far below that item’s marginal.
Splitting the sample into subgroups widens these basins, because the
pile-up makes up a larger share of each smaller sample, so a model that
looks clean on the whole sample deserves the check again on the parts.

## Missing data needs nothing from you

Any indicator containing `NA` is estimated by full-information maximum
likelihood under the usual missing-at-random assumption, selected
automatically. Cases missing on *every* indicator are dropped, with a
warning saying how many, because they contribute nothing while still
counting towards the sample size in BIC.

## Multiple groups: what the log-likelihood is over

With `group =`, mixtureEM maximises the likelihood of the indicators
**given** the group. The grouping variable’s own distribution is not
modelled and its proportions are not counted as parameters.

That is a convention, and programs differ. Software that treats a fully
observed grouping variable as a latent class variable instead adds the
group’s multinomial term to the log-likelihood and its `G - 1`
proportions to the parameter count. Both differences are fixed
constants: they cancel in any likelihood-ratio test and in any
comparison among models fitted here. For reading numbers off such output
directly, a `group` fit also carries `metrics$ll_knownclass` and
`metrics$n_params_knownclass`.

## References

Biernacki, C., Celeux, G., & Govaert, G. (2003). Choosing starting
values for the EM algorithm for getting the highest likelihood in
multivariate Gaussian mixture models. *Computational Statistics & Data
Analysis*, *41*(3-4), 561-575.

Chen, J., Tan, X., & Zhang, R. (2008). Inference for normal mixtures in
mean and variance. *Statistica Sinica*, *18*(2), 443-465.

Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their
Application* (ch. 4). Cambridge University Press.

Day, N. E. (1969). Estimating the components of a mixture of normal
distributions. *Biometrika*, *56*(3), 463-474.

Dziak, J. J., Lanza, S. T., & Tan, X. (2014). Effect size, statistical
power and sample size requirements for the bootstrap likelihood ratio
test in latent class analysis. *Structural Equation Modeling*, *21*(4),
534-552.

Galindo Garre, F., & Vermunt, J. K. (2006). Avoiding boundary estimates
in latent class analysis by Bayesian posterior mode estimation.
*Behaviormetrika*, *33*(1), 43-59.

Hathaway, R. J. (1985). A constrained formulation of maximum-likelihood
estimation for normal mixture distributions. *The Annals of Statistics*,
*13*(2), 795-800.

Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the estimation of
growth mixture models. *Psychological Methods*, *11*(1), 36-53.

Kiefer, J., & Wolfowitz, J. (1956). Consistency of the maximum
likelihood estimator in the presence of infinitely many incidental
parameters. *The Annals of Mathematical Statistics*, *27*(4), 887-906.

Lee, T. K., Wickrama, K. A. S., & O’Neal, C. W. (2023). An introduction
to growth mixture models (GMM). In *International Encyclopedia of
Education* (4th ed., Vol. 14, pp. 646-655). Elsevier.

McLachlan, G. J., & Peel, D. (2000). *Finite Mixture Models*. Wiley.

Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions
about latent class analysis. *Translational Issues in Psychological
Science*, *4*(4), 440-461.

Shireman, E., Steinley, D., & Brusco, M. J. (2017). Examining the effect
of initialization strategies on the performance of Gaussian mixture
modeling. *Behavior Research Methods*, *49*(1), 282-293.

Tekle, F. B., Gudicha, D. W., & Vermunt, J. K. (2016). Power analysis
for the bootstrap likelihood ratio test for the number of classes in
latent class models. *Advances in Data Analysis and Classification*,
*10*(2), 209-224.
