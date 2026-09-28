# Link separately fitted mixture models by the three-step

Relates the latent classes of two or more mixture models fitted on the
same cases – for example a latent class analysis of early attitudes and
a growth mixture model of later achievement – with the bias-adjusted
three-step estimator. Each fitted model is step 1. Step 2 assigns every
case to its most likely class under each model and tabulates how often
that assignment misses the true class. Step 3 fits a multinomial
logistic regression of each occasion's latent class on the previous
occasion's, holding those error tables fixed, so the classification
error does not attenuate the association (Vermunt, 2010; Bakk, Tekle &
Vermunt, 2013; Asparouhov & Muthén, 2014).

## Usage

``` r
link_models(
  models,
  predictors_initial = NULL,
  predictors_transition = NULL,
  correction = c("ML", "none"),
  weights = NULL,
  n_init = 10,
  random_state = NULL,
  standard_errors = TRUE
)
```

## Arguments

- models:

  A list of two or more fitted models (from
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
  or similar), in occasion order, each fitted on the same cases in the
  same row order. Names, if given, label the occasions.

- predictors_initial, predictors_transition:

  Optional covariates for the first occasion's class and for the later
  occasions' classes: a numeric vector, matrix or data frame with one
  row per case.

- correction:

  `"ML"` (default) holds the classification error fixed in step 3;
  `"none"` treats the assigned classes as the true ones, the naive
  classify-analyse baseline.

- weights:

  Optional case weights, treated as frequencies.

- n_init:

  Number of starting points for step 3; the first is all coefficients
  zero, the rest random.

- random_state:

  Optional seed for the random starts.

- standard_errors:

  Logical; compute standard errors from the observed information of
  step 3. They treat the error tables as known, so they do not include
  the uncertainty of the step-1 estimates.

## Value

An object of class `linked_model` with the step-3 log-likelihood
(`loglik`), the number of parameters, BIC, the coefficients and their
covariance matrix, the implied initial class proportions (`initial`) and
transition probabilities (`transitions`, one K_t x K\_(t+1) matrix per
pair of occasions, averaged over cases), the classification-error tables
(rows the assigned class, columns the true class) and the modal
assignments.

## Details

The models may differ in family and in number of classes. Occasion 1's
class is regressed on `predictors_initial`; each later occasion's class
on the previous occasion's class and on `predictors_transition`, with
one slope per destination class shared by every origin class. The last
class is the reference throughout.

## References

Vermunt, J. K. (2010). Latent class modeling with covariates: Two
improved three-step approaches. *Political Analysis*, 18(4), 450-469.

Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
association between latent class membership and external variables using
bias-adjusted three-step approaches. *Sociological Methodology*, 43(1),
272-311.

Asparouhov, T., & Muthén, B. (2014). Auxiliary variables in mixture
modeling: Three-step approaches using Mplus. *Structural Equation
Modeling*, 21(3), 329-341.

Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
transition mixture model using the three-step specification. *Structural
Equation Modeling*, 21(3), 439-454.

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 400
s1 <- sample(1:3, n, replace = TRUE)
s2 <- ifelse(runif(n) < c(0.8, 0.5, 0.2)[s1], 1, 2)
X1 <- sapply(1:5, function(j) rbinom(n, 1, c(0.9, 0.5, 0.1)[s1]))
X2 <- sapply(1:4, function(j) rnorm(n, c(0, 3)[s2]))
a <- fit_mixture(X1, n_classes = 3, measurement = "binary", n_init = 5)
b <- fit_mixture(X2, n_classes = 2, measurement = "continuous", n_init = 5)
link <- link_models(list(early = a, late = b))
link
#> Three-step link of 2 mixture models (early, late); classes 3 -> 2
#> Log-likelihood -677.560, 5 parameters, BIC 1385.077
#> 
#>                              Estimate    SE
#> early Class 1 ~ (Intercept)     0.230 0.145
#> early Class 2 ~ (Intercept)     0.056 0.226
#> late Class 1 ~ (Intercept)     -1.389 0.351
#> late Class 1 ~ early Class 1    2.519 0.401
#> late Class 1 ~ early Class 2    1.371 0.540
#> 
#> Transition probabilities (averaged over cases):
#> early -> late 
#>               late Class 1 late Class 2
#> early Class 1        0.756        0.244
#> early Class 2        0.495        0.505
#> early Class 3        0.200        0.800
#> 
#> Standard errors treat the classification-error tables as known.
# }
```
