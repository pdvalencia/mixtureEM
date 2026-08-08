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

Read that line as a stability check. If the best solution was found by
most of the starts, the surface has been mapped and more starts will not
change the answer. If it was found **once**, raise `n_init` and refit: a
maximum seen a single time may just be the best of a small sample of the
surface, and a different seed could beat it. This is the standard
multi-start report in the latent class literature (Hipp & Bauer, 2006;
Nylund-Gibson & Choi, 2018).

Set `random_state` to make the search reproducible.

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
  is `variances = n_classes` — doubling it if the warning persists. The
  warning prints that number for your model.

Then check the result rather than just that the warning stopped: the
flagged variance should no longer be far below the others, and its class
mean should have come off the floor or ceiling of the response scale.

Theory settles the *form* of this prior and not its size. Any penalty
that diverges as a variance approaches zero, and grows more slowly than
the sample, gives a consistent estimator whatever its constant (Chen,
Tan & Zhang, 2008). The choice of constant is a finite-sample judgement,
which is why the default is deliberately weak and the advice above is a
rule with a check rather than a number to copy.

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

Chen, J., Tan, X., & Zhang, R. (2008). Inference for normal mixtures in
mean and variance. *Statistica Sinica*, *18*(2), 443-465.

Day, N. E. (1969). Estimating the components of a mixture of normal
distributions. *Biometrika*, *56*(3), 463-474.

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

McLachlan, G. J., & Peel, D. (2000). *Finite Mixture Models*. Wiley.

Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions
about latent class analysis. *Translational Issues in Psychological
Science*, *4*(4), 440-461.
