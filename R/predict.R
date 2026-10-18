# ==============================================================================
# Out-of-sample scoring: predict() for new cases
# ==============================================================================

#' Classify new cases from a fitted model
#'
#' Computes the posterior class probabilities of cases that were not in the
#' estimation sample, holding every parameter of the fitted model fixed. This is
#' the E-step of the fitted model applied to new rows: by Bayes' rule a case's
#' posterior is its class-specific likelihood times its prior class probability,
#' normalised over classes. Vermunt and Magidson (2025) show that for the
#' standard emissions this is exactly a logistic (or, with class-specific
#' variances, quadratic) scoring equation, and that missing items change only
#' its constant; running the E-step gives the same numbers without assembling
#' the equation, missing values included.
#'
#' @details
#' **Matching the new data to the fit.** Indicators are matched to the fitted
#' model's columns by name, so `newdata` may hold them in any order and may
#' carry other columns. A factor or character indicator is coded against the
#' levels it had when the model was fitted, never against the levels in
#' `newdata`, and a value the fitted model never saw is an error rather than a
#' silent mis-coding. A binary indicator that was recoded to 0/1 at fit time is
#' recoded the same way. The same holds for a factor or character predictor:
#' a level the model was not fitted with is an error, where scoring it would
#' silently treat the case as belonging to the reference level.
#'
#' **Missing values.** A missing answer simply drops out of that case's
#' likelihood (full-information scoring), whether or not the estimation sample
#' had any. A case missing every indicator has no information to classify on
#' and gets `NA`, the same rule that dropped such rows from the estimation.
#'
#' **Predictors and distal outcomes.** For a model with `predictors` (one step
#' or three steps), the prior class probabilities are the fitted multinomial
#' regression evaluated at the new case's predictors, so `newdata` must carry
#' those columns; a factor predictor is supplied as the factor itself. A case
#' with a missing predictor gets `NA`. A distal outcome is never used: a new
#' case is classified from its indicators only, which is what keeps the outcome
#' out of the classification (Vermunt and Magidson, 2025, sec. 1.4). For a
#' model fitted with a distal outcome in one step, `predict()` on the
#' estimation data therefore differs from [`class_assignments()`], which
#' conditions on the outcome too. It differs for a three-step model with
#' predictors as well: there [`class_assignments()`] reports the step-3
#' posterior, which conditions on the step-1 class assignment, while
#' `predict()` conditions on the indicators themselves, since a new case has
#' no step-1 assignment. On a one-step model the two agree exactly.
#'
#' **Not covered.** Latent transition models (scoring a sequence of statuses
#' needs a choice between filtered and smoothed probabilities), growth models
#' fitted from long data (`fit_lcga()`, `fit_gmm()`, `fit_rmlca()`),
#' multiple-group models, models with direct effects of covariates on items
#' (`predictors_items`), and models combining predictors with a distal outcome
#' in one structural model are refused with an error.
#'
#' @param object A fitted model from [`fit_mixture()`].
#' @param newdata A data frame or matrix holding the indicators (and the
#'   predictors, for a model that has them) of the cases to classify.
#' @param type What to return, as in [`class_assignments()`]: `"modal"` the
#'   most probable class, `"posterior"` the matrix of posterior probabilities,
#'   `"both"` a data frame of the two.
#' @param ... Unused.
#'
#' @return One row (or element) per row of `newdata`, in the format
#'   [`class_assignments()`] returns; `NA` for a case that cannot be scored.
#'
#' @references
#' Vermunt, J. K., & Magidson, J. (2025). Linear logistic scoring equations for
#' latent class and latent profile models: A simple method for classifying new
#' cases. \emph{Structural Equation Modeling}, \emph{32}(3), 541–549.
#' \doi{10.1080/10705511.2024.2375736}
#'
#' @seealso [`class_assignments()`] for the cases in the estimation sample.
#'
#' @examples
#' set.seed(1)
#' X   <- matrix(rbinom(600, 1, rep(c(0.2, 0.8), each = 300)), ncol = 6,
#'               dimnames = list(NULL, paste0("u", 1:6)))
#' fit <- fit_mixture(X[1:80, ], n_classes = 2, measurement = "binary")
#' new <- X[81:100, ]
#' new[1, 2] <- NA                       # a missing answer is fine
#' predict(fit, new)
#' head(predict(fit, new, type = "both"))
#'
#' @export
predict.mixture_model <- function(object, newdata,
                                  type = c("modal", "posterior", "both"), ...) {
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata))
    stop("`newdata` is required; for the estimation sample use ",
         "class_assignments().", call. = FALSE)
  .predict_refuse(object)

  X  <- .predict_indicators(object, newdata)
  ms <- object
  if (anyNA(X)) ms$mm <- .with_missing_emission(ms$mm)

  # Only a class-membership regression carries over to a new case; a distal
  # outcome is not available for one and must not shape its class.
  Y <- NULL
  if (inherits(ms$sm, "covariate")) Y <- .predict_covariates(object, newdata)
  else ms$sm <- NULL

  ok <- rowSums(!is.na(X)) > 0L
  if (!is.null(Y)) {
    bad <- ok & !stats::complete.cases(Y)
    if (any(bad))
      warning(sprintf(paste0("%d case%s with a missing predictor cannot be ",
                             "scored and %s NA."), sum(bad),
                      if (sum(bad) == 1L) "" else "s",
                      if (sum(bad) == 1L) "gets" else "get"), call. = FALSE)
    ok <- ok & !bad
  }

  K    <- object$n_components
  resp <- matrix(NA_real_, nrow(X), K,
                 dimnames = list(NULL, paste("Class", seq_len(K))))
  if (any(ok)) {
    lr <- e_step(ms, X[ok, , drop = FALSE],
                 if (is.null(Y)) NULL else Y[ok, , drop = FALSE])$log_resp
    resp[ok, ] <- exp(lr)
  }

  modal <- rep(NA_integer_, nrow(X))
  modal[ok] <- max.col(resp[ok, , drop = FALSE], ties.method = "first")
  if (type == "modal") return(modal)
  if (type == "posterior") return(resp)
  data.frame(class       = modal,
             probability = resp[cbind(seq_len(nrow(resp)), modal)],
             resp,
             check.names = FALSE)
}

#' @export
predict.lta_model <- function(object, ...)
  stop("predict() does not cover latent transition models yet: scoring a new ",
       "case's sequence of statuses needs a choice between filtered and ",
       "smoothed probabilities. Use class_assignments() for the estimation ",
       "sample.", call. = FALSE)

# The fits whose new-data scoring is not built, each with what is missing.
.predict_refuse <- function(object, what = "predict()") {
  if (inherits(object, c("lcga", "gmm", "rmlca")))
    stop(what, " does not cover growth models fitted from long data yet.",
         call. = FALSE)
  if (!is.null(object$group_info) || inherits(object$mm, "group_blocks") ||
      inherits(object$sm, "group_prevalence"))
    stop(what, " does not cover multiple-group models yet.", call. = FALSE)
  if (.has_dif(object$mm))
    stop(what, " does not cover models with direct effects of covariates on ",
         "items (`predictors_items`).", call. = FALSE)
  if (inherits(object$sm, "nested"))
    stop(what, " does not cover a structural model that combines predictors ",
         "with a distal outcome; fit the outcome with add_outcome() instead.",
         call. = FALSE)
  if (is.null(object$data))
    stop("`object` carries no estimation data to match new cases against.",
         call. = FALSE)
  invisible(TRUE)
}

.has_dif <- function(mm) {
  if (!is.null(mm$models)) return(any(vapply(mm$models, .has_dif, logical(1))))
  any(grepl("_dif$", class(mm)))
}

# A model fitted to complete data has the complete-data likelihood, which
# returns NA on a missing cell. Its missing-data counterpart (`.nan_variant`)
# reads the same parameters and drops the cell instead.
.with_missing_emission <- function(mm) {
  if (!is.null(mm$models)) {
    mm$models <- lapply(mm$models, .with_missing_emission)
    return(mm)
  }
  cl <- class(mm)[1L]
  if (cl %in% names(.nan_variant)) class(mm)[1L] <- .nan_variant[[cl]]
  mm
}

# The new cases' indicators, coded the way the fitted model's columns were.
.predict_indicators <- function(object, newdata) {
  nm <- colnames(object$data)
  if (is.null(nm))
    stop("The fitted model has no indicator names to match `newdata` ",
         "against; refit it with named indicator columns.", call. = FALSE)
  nd <- as.data.frame(newdata, stringsAsFactors = FALSE)
  absent <- setdiff(nm, names(nd))
  if (length(absent))
    stop("`newdata` lacks indicator column(s) the model was fitted on: ",
         paste(absent, collapse = ", "), ".", call. = FALSE)

  levs  <- object$indicator_levels
  X     <- matrix(NA_real_, nrow(nd), length(nm), dimnames = list(NULL, nm))
  for (j in nm) {
    v   <- nd[[j]]
    lev <- levs[[j]]
    if (is.null(lev)) {
      if (!is.numeric(v) && !all(is.na(v)))
        stop(sprintf("Indicator '%s' was numeric when the model was fitted ",
                     j), "and must be numeric in `newdata`.", call. = FALSE)
      X[, j] <- as.numeric(v)
    } else if (is.logical(v)) {
      X[, j] <- as.numeric(v)            # data.matrix() codes FALSE/TRUE as 0/1
    } else {
      if (is.numeric(v) && !all(is.na(v)))
        stop(sprintf(paste0("Indicator '%s' was a factor when the model was ",
                            "fitted; supply it as one in `newdata` so its ",
                            "levels can be matched."), j), call. = FALSE)
      code <- match(as.character(v), lev)
      unseen <- unique(as.character(v)[is.na(code) & !is.na(v)])
      if (length(unseen))
        stop(sprintf("Indicator '%s' has value(s) the model never saw: %s.",
                     j, paste(unseen, collapse = ", ")), call. = FALSE)
      X[, j] <- code
    }
  }

  for (r in object$binary_recode) {
    if (is.null(r$codes))
      stop("This model was fitted by an older version of the package, which ",
           "did not record its binary recoding; refit it to use predict().",
           call. = FALSE)
    v <- X[, r$item]
    if (any(!is.na(v) & !v %in% r$codes))
      stop(sprintf("Binary indicator '%s' has values other than its two.",
                   r$item), call. = FALSE)
    X[, r$item] <- ifelse(is.na(v), NA_real_, as.numeric(v == r$codes[2L]))
  }
  X
}

# The new cases' predictors in the fitted regression's columns. A factor
# predictor was dummy-coded at fit time as `<name>.<level>`, which
# `covariate_terms` lets us rebuild from the factor itself.
.predict_covariates <- function(object, newdata) {
  cn    <- colnames(object$sm$parameters$beta)[-1L]
  terms <- attr(object$Y, "covariate_terms") %||% cn
  nd    <- as.data.frame(newdata, stringsAsFactors = FALSE)
  Y     <- matrix(NA_real_, nrow(nd), length(cn), dimnames = list(NULL, cn))
  for (i in seq_along(cn)) {
    tm <- terms[i]
    if (!tm %in% names(nd))
      stop(sprintf("`newdata` lacks predictor '%s'.", tm), call. = FALSE)
    v <- nd[[tm]]
    if (cn[i] == tm) {
      if (!is.numeric(v))
        stop(sprintf("Predictor '%s' must be numeric in `newdata`.", tm),
             call. = FALSE)
      Y[, i] <- v
    } else {
      # All dummies 0 is the reference level, so a level the fit never saw
      # would silently be scored as the reference; refuse it instead.
      seen   <- attr(object$Y, "covariate_levels")[[tm]]
      unseen <- setdiff(unique(as.character(v[!is.na(v)])), seen)
      if (!is.null(seen) && length(unseen))
        stop(sprintf(paste0("Predictor '%s' has level%s %s in `newdata` ",
                            "that the model was not fitted with."), tm,
                     if (length(unseen) == 1L) "" else "s",
                     paste0("'", unseen, "'", collapse = ", ")),
             call. = FALSE)
      level  <- substring(cn[i], nchar(tm) + 2L)
      Y[, i] <- ifelse(is.na(v), NA_real_, as.numeric(as.character(v) == level))
    }
  }
  Y
}

# ==============================================================================
# Simulating data from a fitted model
# ==============================================================================

#' Simulate data sets from a fitted model
#'
#' Draws new data sets from a fitted model's measurement model, taking its
#' estimates as the truth: each case's class is drawn from the class
#' probabilities, then its indicators from that class's distribution. Refitting
#' the model to many such data sets and counting how often the classes and
#' parameters are recovered is the Monte Carlo route to a sample-size or power
#' question, with the fitted (or a hypothesised) solution as the population.
#'
#' @details
#' Classes are drawn from the fitted class proportions. For a model with
#' predictors they are drawn from each case's fitted class-membership
#' probabilities given its observed predictors, which are held fixed rather
#' than resampled, so `n` must equal the estimation sample size; the same rule
#' holds for a growth mixture model with covariates on the growth factors. A
#' distal outcome is not simulated. The data come back complete and on the
#' model's own coding (0/1 for a binary item, 1 to M for a categorical one),
#' which [`fit_mixture()`] reads back unchanged.
#'
#' To simulate from a hypothesised solution rather than a fitted one, fit any
#' model of the right shape and overwrite its parameters before calling
#' `simulate()`. Latent transition models, multiple-group models and models
#' with direct effects on items are not covered.
#'
#' @param object A fitted model from [`fit_mixture()`], [`fit_lcga()`] or
#'   [`fit_gmm()`].
#' @param nsim Number of data sets.
#' @param seed `NULL`, or an integer passed to [set.seed()]; the random number
#'   state is restored afterwards, as for [stats::simulate()].
#' @param n Cases per data set; defaults to the estimation sample size.
#' @param ... Unused.
#'
#' @return A list of `nsim` data frames named like the fitted indicators, each
#'   carrying the drawn classes as its `"classes"` attribute, with the random
#'   number state as the list's `"seed"` attribute.
#'
#' @seealso [`predict.mixture_model()`] to classify the simulated cases.
#'
#' @examples
#' set.seed(1)
#' X   <- matrix(rbinom(600, 1, rep(c(0.2, 0.8), each = 300)), ncol = 6,
#'               dimnames = list(NULL, paste0("u", 1:6)))
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#' sims <- simulate(fit, nsim = 2, seed = 42, n = 200)
#' head(sims[[1]])
#' table(attr(sims[[1]], "classes"))
#'
#' @importFrom stats simulate
#' @export
simulate.mixture_model <- function(object, nsim = 1, seed = NULL, n = NULL,
                                   ...) {
  .simulate_refuse(object)
  N <- nrow(object$data)
  n <- n %||% N
  covariate <- inherits(object$sm, "covariate")
  if (covariate && n != N)
    stop(sprintf(paste0("A model with predictors simulates classes ",
                        "conditional on the observed predictors, so a data ",
                        "set must have one case per original case (%d ",
                        "requested, %d available)."), n, N), call. = FALSE)

  # The seed contract of stats::simulate(): a supplied seed is used and the
  # caller's random number state put back afterwards.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    stats::runif(1)
  if (is.null(seed)) {
    rng <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    saved <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv))
    set.seed(seed)
    rng <- structure(seed, kind = as.list(RNGkind()))
  }

  K     <- object$n_components
  probs <- if (covariate) exp(log_likelihood(object$sm, object$Y)) else NULL
  nm    <- colnames(object$data)
  out <- lapply(seq_len(nsim), function(s) {
    cls <- if (covariate)
      vapply(seq_len(n), function(i) sample.int(K, 1L, prob = probs[i, ]),
             integer(1))
    else sample.int(K, n, replace = TRUE, prob = object$weights)
    d <- as.data.frame(generate_synthetic_data(object$mm, cls, n))
    if (!is.null(nm)) names(d) <- nm
    # "class" is already the data frame's own attribute, so the drawn classes
    # ride as "classes".
    attr(d, "classes") <- cls
    d
  })
  attr(out, "seed") <- rng
  out
}

#' @export
simulate.lta_model <- function(object, ...)
  stop("simulate() does not cover latent transition models yet.",
       call. = FALSE)

.simulate_refuse <- function(object) {
  if (!is.null(object$group_info) || inherits(object$mm, "group_blocks") ||
      inherits(object$sm, "group_prevalence"))
    stop("simulate() does not cover multiple-group models yet.", call. = FALSE)
  if (.has_dif(object$mm))
    stop("simulate() does not cover models with direct effects of ",
         "covariates on items (`predictors_items`).", call. = FALSE)
  if (is.null(object$data))
    stop("`object` carries no estimation data to size the simulation by.",
         call. = FALSE)
  invisible(TRUE)
}
