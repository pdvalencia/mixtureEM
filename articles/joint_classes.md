# Two Latent Class Variables: Association, Covariates and Moderation

``` r

library(mixtureEM)
```

Some research questions are about two typologies at once. Do the
profiles of childhood adversity a survey finds go with the profiles of
substance use it also finds? Is that link the same for girls and boys?
Each typology is a latent class variable measured by its own items, and
the question is about the *association* between the two.
[`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md)
fits both latent variables in one model, with the association written as
a multinomial logistic regression of one on the other, the loglinear
model with latent variables of Hagenaars (1993). Covariates can predict
either variable, and a covariate’s effect on the second can differ
across the classes of the first, which is how a moderated association is
written.

What a latent class analysis is, and how to choose the number of
classes, is in
[`vignette("class_enumeration")`](https://pdvalencia.github.io/mixtureEM/articles/class_enumeration.md);
here each variable’s number of classes is taken as settled.

## The data

The example is simulated, so that the answer is known. For 1,500
adolescents there are six yes/no items about adverse childhood
experiences (`a1` to `a6`) and five items about substance use, each
answered never (1), sometimes (2) or often (3) (`s1` to `s5`). Behind
them are two latent variables with three classes each: adversity (`ca`:
low, household adversity, high) and substance use (`cs`: abstains,
drinks, several substances). High adversity goes with several
substances, and being a girl lowers the odds of several substances, most
of all after high adversity. Three answers in a hundred are missing at
random.

``` r

set.seed(7070)
n <- 1500
female <- rbinom(n, 1, 0.5)
# Adversity class: 1 low, 2 household, 3 high; girls a little more often 2.
lp <- cbind(0, -0.6 + 0.5 * female, -1.2)
ca <- apply(exp(lp) / rowSums(exp(lp)), 1, function(p) sample(1:3, 1, prob = p))
# Substance-use class given adversity: 1 abstains, 2 drinks, 3 several
# substances. Logits against class 1, by adversity class; being female lowers
# the odds of several substances, most of all after high adversity.
assoc <- rbind(c(0, -0.4, -2.2), c(0, 0.3, -0.6), c(0, 0.6, 0.9))
lc <- assoc[ca, ] + cbind(0, 0, -c(0.2, 0.4, 1.2)[ca] * female)
cs <- apply(exp(lc) / rowSums(exp(lc)), 1, function(p) sample(1:3, 1, prob = p))
# Six yes/no adversity items and five never/sometimes/often use items.
p_a <- rbind(c(.10, .08, .12, .05, .10, .06),
             c(.70, .65, .55, .20, .15, .10),
             c(.80, .75, .70, .65, .60, .55))
adversity <- matrix(rbinom(n * 6, 1, p_a[ca, ]), n, 6,
                    dimnames = list(NULL, paste0("a", 1:6)))
use_level <- rbind(c(0, 0, 0, 0, 0), c(2, 0, 1, 0, 0), c(2, 2, 2, 1, 1))
use <- matrix(NA_integer_, n, 5, dimnames = list(NULL, paste0("s", 1:5)))
for (j in 1:5) {
  hi <- use_level[cs, j]
  r <- runif(n)
  use[, j] <- 1 + ifelse(hi == 0, (r > .80) + (r > .95),
              ifelse(hi == 1, (r > .30) + (r > .80), (r > .10) + (r > .40)))
}
adversity[runif(n * 6) < 0.03] <- NA
use[runif(n * 5) < 0.03] <- NA
```

`ca` and `cs` are the truth the analysis will try to recover; a real
analysis never has them. All it has is `adversity`, `use` and `female`.

## One step or two

The model can be estimated in one step, every parameter at once, or in
two (Bakk & Kuha, 2018). In two steps, step 1 fits each latent
variable’s measurement model on its own, as two separate latent class
analyses, and step 2 holds those item parameters fixed while it
estimates the association and the covariate effects. The two-step
estimator has a property that matters when several structural models are
compared: the classes are defined by their own items only, so they mean
the same thing in every model, and a covariate added in step 2 cannot
change what a class is. Its standard errors include the uncertainty of
step 1. The analysis below is two-step, with the one-step fit as a check
at the end.

Every fit switches off the package’s default priors, which keep
estimates away from the boundary and are the right default for
exploratory work. Without them the fits are plain maximum likelihood,
the estimator the two-step method and the likelihood-ratio tests below
are defined for;
[`vignette("estimation")`](https://pdvalencia.github.io/mixtureEM/articles/estimation.md)
explains what the priors add and when to keep them.

``` r

ml <- list(smoothing = 0, bayes_constants = list(categorical = 0))
```

## The association with covariates

The first structural model has the association and `female` predicting
both variables, its effect on `cs` the same in every adversity class
(`cs_effects = "common"`, the default).

``` r

fem <- data.frame(female = female)
m1 <- fit_joint_classes(adversity, use, n_ca = 3, n_cs = 3,
                        measurement_cs = "ordinal",
                        predictors_ca = fem, predictors_cs = fem,
                        n_steps = 2, n_init = 20, random_state = 1,
                        smoothing = ml$smoothing,
                        bayes_constants = ml$bayes_constants)
```

Before the association, what the classes are.
[`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
gives each variable’s item profiles under its own classes, beside the
observed share of each answer:

``` r

measurement_summary(m1)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> ca: item response probabilities by ca class
#> Indicator            | Overall |    CA#1 |    CA#2 |    CA#3
#> ------------------------------------------------------------
#> a1                   |   0.399 |   0.079 |   0.813 |   0.731
#> a2                   |   0.384 |   0.103 |   0.629 |   0.811
#> a3                   |   0.335 |   0.118 |   0.564 |   0.623
#> a4                   |   0.201 |   0.056 |   0.240 |   0.528
#> a5                   |   0.196 |   0.070 |   0.219 |   0.493
#> a6                   |   0.143 |   0.060 |   0.000 |   0.528
#> 
#> cs: item response probabilities by cs class
#> Indicator            | Overall |    CS#1 |    CS#2 |    CS#3
#> ------------------------------------------------------------
#> s1 (Cat 1)           |   0.385 |   0.802 |   0.000 |   0.156
#> s1 (Cat 2)           |   0.260 |   0.198 |   0.313 |   0.306
#> s1 (Cat 3)           |   0.356 |   0.000 |   0.687 |   0.538
#> s2 (Cat 1)           |   0.701 |   0.788 |   0.780 |   0.000
#> s2 (Cat 2)           |   0.161 |   0.164 |   0.143 |   0.221
#> s2 (Cat 3)           |   0.138 |   0.048 |   0.077 |   0.779
#> s3 (Cat 1)           |   0.478 |   0.731 |   0.306 |   0.099
#> s3 (Cat 2)           |   0.337 |   0.203 |   0.496 |   0.261
#> s3 (Cat 3)           |   0.185 |   0.066 |   0.198 |   0.641
#> s4 (Cat 1)           |   0.731 |   0.797 |   0.791 |   0.207
#> s4 (Cat 2)           |   0.198 |   0.146 |   0.174 |   0.519
#> s4 (Cat 3)           |   0.071 |   0.057 |   0.035 |   0.274
#> s5 (Cat 1)           |   0.743 |   0.814 |   0.806 |   0.188
#> s5 (Cat 2)           |   0.191 |   0.139 |   0.144 |   0.600
#> s5 (Cat 3)           |   0.066 |   0.046 |   0.051 |   0.213
#> 
#> A two-category item shows its second category: for an item coded 0/1, the probability of a 1.
#> 
#> At the boundary: a6 in class 2 (ca); s1 in class 2 (cs); s1 in class 1 (cs); s2 in class 3 (cs). These probabilities have run to 0 or 1, so the class is defined partly by an item every case in it gives the same answer to, and their standard errors are not interpretable. There are two ways on: read it substantively, since an item every member of a class answers identically is often the finding rather than a fault; or, if that parameter needs a standard error, refit with a stronger prior than the default of 1 - `bayes_constants = list(categorical = 2)` - which holds the estimate off the edge at the cost of shrinking it slightly toward the item's marginal.
#> =========================================================
```

The adversity classes are the three that were simulated: `CA#1` reports
almost none of the six experiences, `CA#2` the first three (household
adversity) and few of the rest, `CA#3` all six at around one half or
more. The substance-use classes are too: `CS#1` mostly answers never,
`CS#2` drinks often (`s1`) and sometimes uses the third substance, and
`CS#3` uses the first three substances often and the last two sometimes.

Four probabilities have run to 0 or 1, and the note under the tables
says how to read them. The household class never reporting `a6` is a
finding about that class rather than a fault of the fit; it was
simulated at 0.10, and in a class of a few hundred the estimate of a
rare answer can be zero.

Now the association:

``` r

m1
#> Two latent class variables: ca (3 classes, 6 items), cs (3 classes, 5 items)
#> Log-likelihood -10499.379, 60 free parameters, BIC 21437.551
#> 
#> Two-step: item parameters fixed at the separate latent class analyses
#> 
#> Association: log-odds of cs against its last class
#>    cs      term   ca estimate    se      z        p
#>  CS#1 Intercept       -0.0671 0.263 -0.255 7.98e-01
#>  CS#1      CA#1 CA#1   2.3990 0.344  6.964 3.31e-12
#>  CS#1      CA#2 CA#2   0.6601 0.404  1.634 1.02e-01
#>  CS#1    female        0.6574 0.221  2.979 2.89e-03
#>  CS#2 Intercept        0.6824 0.235  2.910 3.62e-03
#>  CS#2      CA#1 CA#1   1.4005 0.329  4.253 2.11e-05
#>  CS#2      CA#2 CA#2   0.0644 0.391  0.165 8.69e-01
#>  CS#2    female        0.5165 0.220  2.349 1.88e-02
#> 
#> Joint class probabilities P(ca, cs)
#>       CS#1  CS#2  CS#3
#> CA#1 0.297 0.218 0.023
#> CA#2 0.103 0.110 0.039
#> CA#3 0.055 0.109 0.045
#> 
#> Conditional probabilities P(cs | ca)
#>       CS#1  CS#2  CS#3
#> CA#1 0.552 0.405 0.042
#> CA#2 0.408 0.436 0.156
#> CA#3 0.264 0.522 0.215
```

Every log-odds is against the last class of each variable:
`CS#1 ON CA#1` of 2.40 says that low adversity raises the odds of
abstaining, rather than using several substances, by a factor of
`exp(2.40) = 11` compared with high adversity. The conditional
probabilities say the same thing on the probability scale: 4% of the
low-adversity class use several substances, 22% of the high-adversity
class.

## Is there an association at all?

The model with `association = "none"` holds the two variables
independent given the covariates. Fitted on the same step 1
(`step1 = m1`), it is the restricted model of a likelihood-ratio test of
the association, on (3 - 1) x (3 - 1) = 4 degrees of freedom.

``` r

m0 <- fit_joint_classes(adversity, use, n_ca = 3, n_cs = 3,
                        measurement_cs = "ordinal",
                        predictors_ca = fem, predictors_cs = fem,
                        association = "none", n_steps = 2, step1 = m1,
                        smoothing = ml$smoothing,
                        bayes_constants = ml$bayes_constants)
lr_test(m0, m1)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -10539.5218   parameters = 56
#>   Full       : LL =  -10499.3791   parameters = 60
#>   -2 x diff  : 80.2854   df = 4   p = 1.515e-16
#>   The restriction is rejected: the full model fits significantly better.
```

The two typologies are clearly associated.

## Is the association the same for girls and boys?

`cs_effects = c(female = "by_ca")` lets the effect of `female` on
substance use differ by adversity class. Against `m1`, that is a test of
whether sex moderates the association, again on four degrees of freedom:
one more slope per `cs` logit for each adversity class but the last.

``` r

m2 <- fit_joint_classes(adversity, use, n_ca = 3, n_cs = 3,
                        measurement_cs = "ordinal",
                        predictors_ca = fem, predictors_cs = fem,
                        cs_effects = c(female = "by_ca"),
                        n_steps = 2, step1 = m1,
                        smoothing = ml$smoothing,
                        bayes_constants = ml$bayes_constants)
lr_test(m1, m2)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -10499.3791   parameters = 60
#>   Full       : LL =  -10494.5294   parameters = 64
#>   -2 x diff  : 9.6996   df = 4   p = 0.0458
#>   The restriction is rejected: the full model fits significantly better.
```

The moderation is significant at the 5% level, if not by much. The
coefficients of a moderated model are hard to read directly (each
`female | CA#k` row is a difference from the last adversity class), so
the result is best reported as the probabilities of each substance-use
class given each adversity class, for boys and for girls.
[`joint_probabilities()`](https://pdvalencia.github.io/mixtureEM/reference/joint_probabilities.md)
computes them at fixed covariate values, with delta-method standard
errors:

``` r

boys  <- joint_probabilities(m2, at = list(female = 0))
girls <- joint_probabilities(m2, at = list(female = 1))
cbind(boys[, c("ca", "cs")], boys = round(boys$estimate, 3),
      girls = round(girls$estimate, 3))
#>     ca   cs  boys girls
#> 1 CA#1 CS#1 0.525 0.590
#> 2 CA#1 CS#2 0.428 0.373
#> 3 CA#1 CS#3 0.047 0.037
#> 4 CA#2 CS#1 0.344 0.447
#> 5 CA#2 CS#2 0.510 0.396
#> 6 CA#2 CS#3 0.146 0.157
#> 7 CA#3 CS#1 0.283 0.245
#> 8 CA#3 CS#2 0.400 0.662
#> 9 CA#3 CS#3 0.317 0.092
```

After high adversity (`CA#3`), 32% of boys are in the several-substances
class and 9% of girls; after low adversity (`CA#1`) both are under 5%.
That is the moderation that was simulated: the link from adversity to
using several substances is much stronger for boys.

Any other function of the coefficients, with its standard error, is
available from
[`derived_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/derived_estimates.md),
using the parameter names
[`joint_coefficients()`](https://pdvalencia.github.io/mixtureEM/reference/joint_coefficients.md)
documents.

## The one-step fit as a check

The same three models in one step, every parameter estimated at once:

``` r

one <- lapply(list(m0 = "none", m1 = "common", m2 = "by_ca"), function(s)
  fit_joint_classes(adversity, use, n_ca = 3, n_cs = 3,
                    measurement_cs = "ordinal",
                    predictors_ca = fem, predictors_cs = fem,
                    association = if (s == "none") "none" else "free",
                    cs_effects = if (s == "by_ca") c(female = "by_ca") else "common",
                    n_init = 20, random_state = 1, standard_errors = FALSE,
                    smoothing = ml$smoothing,
                    bayes_constants = ml$bayes_constants))
lr_test(one$m0, one$m1)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -10539.1580   parameters = 56
#>   Full       : LL =  -10497.3270   parameters = 60
#>   -2 x diff  : 83.6620   df = 4   p = < 1e-16
#>   The restriction is rejected: the full model fits significantly better.
lr_test(one$m1, one$m2)
#> 
#> Likelihood-ratio test for nested models
#> ---------------------------------------------------------
#>   Restricted : LL =  -10497.3270   parameters = 60
#>   Full       : LL =  -10490.8376   parameters = 64
#>   -2 x diff  : 12.9788   df = 4   p = 0.01138
#>   The restriction is rejected: the full model fits significantly better.
```

The conclusions do not change: a strong association, and a moderation by
sex that is significant at the 5% level. The one-step log-likelihoods
are higher, as they must be, since step 1’s item parameters are free
there; the price is that the classes’ meaning can shift a little from
one structural model to the next.

## Summary

- [`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md)
  fits two latent class variables and their association in one model;
  [`measurement_summary()`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
  shows what each variable’s classes are.
- `association = "none"` and `cs_effects = "by_ca"` give the restricted
  and the extended models of the two natural tests, run by
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md).
- With `n_steps = 2`, fit step 1 once and hand it to every structural
  model with `step1 =`, so that all of them share one measurement model.
- Report a moderated association as conditional probabilities at chosen
  covariate values, with
  [`joint_probabilities()`](https://pdvalencia.github.io/mixtureEM/reference/joint_probabilities.md).

## References

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871-892.

Hagenaars, J. A. (1993). *Loglinear models with latent variables*. Sage.
