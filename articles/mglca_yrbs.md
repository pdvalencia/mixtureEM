# Multiple-Group LCA: Do the Classes Differ by Grade?

``` r

library(mixtureEM)
```

## The question

A latent class model assumes the classes mean the same thing to every
case in the data. When the data span groups — here, four grades of high
school — that assumption is worth checking rather than taking for
granted, and it splits into two separate questions. First, *measurement
invariance*: do the twelve health-risk items mean the same thing to a
ninth-grader as to a twelfth-grader, in the sense that the same class
membership implies the same item-response probabilities regardless of
grade? Second, if invariance holds, *prevalence*: are the risk classes
more common in some grades than others? The second question is only
interpretable once the first is settled — a prevalence comparison
between groups whose items do not function the same way is comparing
different things.

This uses the same `yrbs2005` data as
[`vignette("survey_lca")`](https://pdvalencia.github.io/mixtureEM/articles/survey_lca.md),
but unweighted, matching how the textbook this pipeline follows presents
it; see that vignette for the same items analysed with the sampling
design in play.

## Data

`grade` is a factor with levels 9 through 12, with no missing values in
the bundled sample:

``` r

data(yrbs2005)
items <- yrbs2005[, 3:14]
grade <- yrbs2005$grade
table(grade)
#> grade
#>    9   10   11   12 
#> 3332 3470 3529 3509
```

## The baseline: prevalences free, measurement invariant

``` r

set.seed(1)
m_free <- fit_mixture(items, n_classes = 5, measurement = "binary",
                      group = grade, group_effects = "prevalence",
                      n_init = 20, max_iter = 2000, n_steps = 1)
m_free
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 5
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 74 iterations)
#> Missing Data       : 7186 / 166080 cells (4.3%) in 12 items — FIML (MAR assumption)
#> ---------------------------------------------------------
#>   Log-Likelihood : -48032.87
#>   Parameters     : 76
#>   AIC            : 96217.73
#>   BIC            : 96790.42
#>   SABIC          : 96548.90
#>   Rel. Entropy   : 0.8046
#>   Best solution  : found by 6 of 20 starts
#>   (Comparing with software that counts the grouping variable's own proportions? Use metrics$ll_knownclass = -67215.74 and metrics$n_params_knownclass = 79.)
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 66.06%
#>   Class 2: 13.17%
#>   Class 3: 10.13%
#>   Class 4: 5.91%
#>   Class 5: 4.73%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```

`group_effects = "prevalence"` pools the item-response probabilities
across grades — that pooling *is* the measurement-invariance assumption
— while letting each grade have its own class sizes. A four-group,
five-class fit on this many cases is materially more expensive per start
than the package’s other vignettes; the solution above replicated across
6 of the 20 starts.

## The constrained model and the omnibus test

Fitting the same model with `group_effects = "none"` forces every grade
to share the same class sizes as well, which is the null hypothesis that
prevalence does not differ by grade at all:

``` r

set.seed(1)
m_none <- fit_mixture(items, n_classes = 5, measurement = "binary",
                      group = grade, group_effects = "none",
                      n_init = 5, max_iter = 2000, n_steps = 1)
lr_test(m_none, m_free)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -48345.2591   parameters = 64
#>   Full       : LL =  -48032.8661   parameters = 76
#>   -2 x diff  : 624.7859   df = 12   p = < 1e-16
#>   The restriction is rejected: the full model fits significantly better.
```

[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
takes the restricted model first and the full model second — reversing
the order is an error. With four grades and five classes there are
(4 - 1) x (5 - 1) = 12 free prevalence contrasts between the two models,
and the printed `df` above confirms it. The restriction is rejected:
class prevalence does differ by grade. This reproduces the “all five
latent classes” comparison from the textbook’s prevalence-differences
table.

Because `group = grade` is supplied, both fits above also report a
log-likelihood on the known-class scale, `metrics$ll_knownclass` — that
is the number that lines up with what other programs print for this kind
of model, as opposed to `metrics$ll`, which is on this package’s own
1-step scale.

``` r

m_free$metrics$ll_knownclass
#> [1] -67215.74
m_none$metrics$ll_knownclass
#> [1] -67528.13
```

## Per-class tests: does *this* class’s prevalence differ by grade?

The omnibus test above says prevalence differs by grade somewhere among
the five classes, but the textbook goes further and asks the question
separately for each class: hold Low Risk’s share fixed across grades
while the other four float and renormalise, refit, and see how much
likelihood that costs. Repeated five times, this reproduces the
chapter’s Table 5.24. `group_prevalence_equal` fits exactly this
restriction — pass it the index of the class to freeze and it holds that
column of the per-group class probabilities to one shared value while
the rest adjust freely.

There is a wrinkle a straightforward call does not handle. A mixture
model’s classes are not identified independently across separate fits:
nothing pins “class 1” to the same substantive group of respondents from
one
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
call to the next, so a search run five times, once per requested class,
is free to relabel *which* real class lands in whichever slot was asked
for. Left to its own devices the optimiser takes the escape hatch every
time — it reports whichever class is *cheapest* to freeze, regardless of
which one was requested, because that fit has the higher likelihood and
wins the random-restart comparison. Running `group_prevalence_equal`
directly, five times, mostly reproduces the same answer (the cost of
freezing High Risk, the cheapest class) instead of five different ones.

The fix is to anchor each fit to `m_free`’s own solution rather than
searching from scratch.
[`class_sizes()`](https://pdvalencia.github.io/mixtureEM/reference/class_sizes.md)’s
`by_group` attribute gives the per-grade class probabilities `m_free`
already settled on; seeding a fit at those values, plus `m_free`’s own
item-response probabilities, starts EM in the one basin where “class k”
already means the same thing it means in `m_free`.
[`fit_mixture_internal()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture_internal.md)
(the lower-level engine behind
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md))
takes a `warm_start` function for exactly this, and `n_init = 0` tells
it to refine only from that seed rather than compete it against random
restarts:

``` r

by_group <- attr(class_sizes(m_free), "by_group")
n_grades <- nlevels(grade)
n_classes <- 5
gamma0 <- matrix(0, n_grades, n_classes)
for (i in seq_len(nrow(by_group))) {
  gr <- match(by_group$group[i], levels(grade))
  gamma0[gr, by_group$class[i]] <- by_group$proportion[i]
}

anchor_at_m_free <- function(model_state, X, Y) {
  model_state$weights <- m_free$weights
  model_state$mm <- m_free$mm
  model_state$sm$parameters$gamma <- gamma0
  model_state
}

class_labels <- c("Low Risk", "Binge Drinkers", "Early Experimenters",
                  "High Risk", "Sexual Risk-Takers")
freeze_class <- function(k) {
  fit_mixture_internal(X = as.matrix(items), Y = as.integer(grade),
                       n_components = n_classes, measurement = "binary",
                       structural = "group_prevalence", n_steps = 1,
                       n_init = 0, max_iter = 2000,
                       warm_start = anchor_at_m_free,
                       n_groups = n_grades, frozen = k)
}
```

``` r

per_class <- lapply(1:5, function(k) lr_test(freeze_class(k), m_free))
data.frame(class = class_labels,
          dG2 = vapply(per_class, function(t) t$statistic, numeric(1)),
          df = vapply(per_class, function(t) t$df, numeric(1)),
          p_value = vapply(per_class, function(t) t$p_value, numeric(1)))
#>                 class         dG2 df       p_value
#> 1            Low Risk  61.2464552  3  3.183513e-13
#> 2      Binge Drinkers 468.0697346  3 3.962496e-101
#> 3 Early Experimenters 176.4318330  3  5.199765e-38
#> 4           High Risk   0.4150451  3  9.371173e-01
#> 5  Sexual Risk-Takers   5.6799906  3  1.282609e-01
```

Each row matches the textbook’s Table 5.24 to about one decimal (61.3,
468.5, 176.6, 0.4 and 5.7, in the same class order as above), and the
pattern makes substantive sense against the item plot above: Binge
Drinkers, the class whose share changes most across grades in that plot,
is the most expensive to force flat by a wide margin, while High Risk
barely moves across grades and costs almost nothing to freeze. The
all-classes restriction from the previous section (dG2 = 624.8, df = 12,
against the book’s 625.5) is the sum of a similar idea applied to every
class’s share at once, not to each one separately.

## Testing measurement invariance itself

Section 3’s pooling of the item probabilities across grades was an
assumption, not a given — it can be tested the same way the prevalence
question was. `group_effects = "both"` frees the item-response
probabilities by grade as well as the class sizes, and comparing it with
the invariant baseline via
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
tests exactly that assumption:

``` r

set.seed(1)
m_both <- fit_mixture(items, n_classes = 5, measurement = "binary",
                      group = grade, group_effects = "both",
                      n_init = 20, max_iter = 2000, n_steps = 1)
#> Warning: The reported solution was found by 1 of 21 starts that ran to
#> convergence, out of 20 requested. EM climbs the peak it starts nearest, so a
#> maximum seen once may be the best of a small sample of the likelihood surface
#> rather than the best there is: refit with n_init = 100 before reporting.
lr_test(m_free, m_both)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -48032.8661   parameters = 76
#>   Full       : LL =  -47528.9757   parameters = 256
#>   -2 x diff  : 1007.7808   df = 180   p = < 1e-16
#>   The restriction is rejected: the full model fits significantly better.
```

`m_both` carries a warning worth repeating rather than hiding: only 1 of
the 20 requested starts converged to the log-likelihood reported above,
so this solution has not been shown to be the global optimum the way
`m_free`’s was. 256 parameters spread across five classes and four
grades is a much harder surface to search than the models earlier in
this vignette, and `n_init = 20` is not enough to be confident of it
here — the package’s own advice is to refit with `n_init = 100` before
reporting this specific comparison as a finding rather than an
illustration.

With N in the thousands, a likelihood-ratio test like this one also has
power to detect even a trivially small departure from invariance, so a
rejection on its own does not settle the question — it is worth reading
the information criteria alongside the test rather than the p-value in
isolation. Here the test rejects invariance decisively, which on its own
would be unsurprising at this sample size. But the two information
criteria disagree about which model is preferable: AIC favours the fully
group-varying model (it is smaller for `m_both`), while BIC — which
penalises the 180 extra parameters more heavily — favours the invariant
`m_free`. That split is itself informative: the departures from
invariance are real enough for the likelihood-ratio test to see them at
this N, but not so large that a criterion penalising complexity thinks
they are worth the extra parameters. A researcher would want to look at
*which* items differ by grade — and would want the better-replicated fit
above — before deciding whether to treat the invariant model as good
enough to interpret.

## Labelling

``` r

measurement_summary(m_free)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator             | Overall | Class 1 | Class 2 | Class 3 | Class 4 | Class 5
#> --------------------------------------------------------------------------------- 
#> smoked_before_13      |   0.154 |   0.037 |   0.125 |   0.630 |   0.639 |   0.256
#> smoked_daily_30d      |   0.120 |   0.018 |   0.282 |   0.246 |   0.653 |   0.169
#> drove_drinking        |   0.105 |   0.005 |   0.442 |   0.107 |   0.452 |   0.131
#> first_drink_before_13 |   0.255 |   0.134 |   0.175 |   0.772 |   0.683 |   0.429
#> binge_drink_30d       |   0.247 |   0.078 |   0.744 |   0.438 |   0.788 |   0.212
#> marijuana_before_13   |   0.089 |   0.005 |   0.026 |   0.379 |   0.549 |   0.263
#> cocaine_ever          |   0.083 |   0.004 |   0.181 |   0.069 |   0.838 |   0.030
#> glue_ever             |   0.116 |   0.053 |   0.151 |   0.256 |   0.564 |   0.037
#> meth_ever             |   0.058 |   0.004 |   0.087 |   0.030 |   0.692 |   0.008
#> ecstasy_ever          |   0.061 |   0.004 |   0.100 |   0.048 |   0.624 |   0.066
#> sex_before_13         |   0.074 |   0.016 |   0.001 |   0.141 |   0.288 |   0.704
#> sex_4plus_partners    |   0.170 |   0.060 |   0.312 |   0.125 |   0.558 |   0.927
#> 
#> Missing data: 7186 of 166080 cells (4.3%) across 12 items, handled via FIML (MAR assumption).
#> 
#> At the boundary: sex_before_13 in class 2. These probabilities have run to 0 or 1, so the class is defined partly by an item every case in it gives the same answer to, and their standard errors are not interpretable. There are two ways on: read it substantively, since an item every member of a class answers identically is often the finding rather than a fault; or, if that parameter needs a standard error, refit with a stronger prior than the default of 1 - `bayes_constants = list(categorical = 2)` - which holds the estimate off the edge at the cost of shrinking it slightly toward the item's marginal.
#> =========================================================
```

This fit is unweighted and grouped, and its classes do not read the same
way
[`vignette("survey_lca")`](https://pdvalencia.github.io/mixtureEM/articles/survey_lca.md)’s
do — which is precisely the point of comparing them. That fit had one
class high on essentially every item; this one splits that pattern in
two: class 4 is where the hard-drug items (`cocaine_ever`, `glue_ever`,
`meth_ever`, `ecstasy_ever`) and the smoking items peak, while class 5
is instead the class where the two sexual-risk items (`sex_before_13`,
`sex_4plus_partners`) are far higher than in any other class, at
moderate levels of everything else. There is no single “high risk across
everything” class in this solution. In size order, the five classes are
**Low Risk, Binge Drinkers, Early Experimenters, High Risk,** and
**Sexual Risk-Takers**:

``` r

plot(m_free,
     class_labels = c("Low Risk", "Binge Drinkers", "Early Experimenters",
                      "High Risk", "Sexual Risk-Takers"),
     main = "Health-risk behavior classes by grade (YRBS 2005, unweighted)")
```

![](mglca_yrbs_files/figure-html/plot-1.png)

Class 2 (Binge Drinkers) stands out on `drove_drinking` and
`binge_drink_30d` specifically; class 3 (Early Experimenters) stands out
on the two early-onset items, `smoked_before_13` and
`first_drink_before_13`, more than on the items that come later in
adolescence. High Risk and Sexual Risk-Takers replace the single “High
risk” class from the weighted fit because the hard-drug cluster and the
sexual-risk cluster no longer travel together in one class here.

## References

Collins, L. M., & Lanza, S. T. (2010). *Latent class and latent
transition analysis: With applications in the social, behavioral, and
health sciences*. Wiley.
