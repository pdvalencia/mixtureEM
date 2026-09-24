# ==============================================================================
# add_covariates() / add_outcome(): stepwise analyses on a fitted model
#
# Both verbs take a fitted mixture_model and run only steps 2-3 of the
# three-step approach on its stored step-1 solution. The measurement model is
# never re-estimated, so the classes the user inspected are exactly the classes
# the structural model describes.
# ==============================================================================

# Common validation for both verbs. Returns the fit with any previous
# structural model cleared and measurement-only posteriors restored.
.check_stepwise_fit <- function(fit, verb) {
  if (inherits(fit, "lta_model"))
    stop("`", verb, "()` does not support LTA models. Supply covariates ",
         "through `fit_lta()`'s `predictors_initial` / `predictors_transition` ",
         "arguments instead.", call. = FALSE)
  if (!inherits(fit, "mixture_model"))
    stop("`fit` must be a fitted model returned by fit_mixture() (or ",
         "fit_rmlca(), fit_lcga(), fit_gmm()).", call. = FALSE)
  if (!is.null(fit$sm) && fit$n_steps == 1L)
    stop("This model was fit in one step, so its classes already condition ",
         "on the structural model and cannot be reused as a step-1 solution. ",
         "Fit the unconditional model first (fit_mixture() without ",
         "`predictors` or `outcome`), then call `", verb, "()` on it.",
         call. = FALSE)
  if (!is.null(fit$group_effects) &&
      fit$group_effects %in% c("both", "prevalence"))
    stop("This model already uses `group` as a class-membership predictor ",
         "(group_effects = \"", fit$group_effects, "\"). Combine the group ",
         "and the new predictors in a single fit_mixture(group = , ",
         "predictors = ) call instead.", call. = FALSE)

  if (!is.null(fit$sm)) {
    message("Replacing the existing structural model; the step-1 measurement ",
            "model is reused unchanged.")
    # The ML correction overwrites $log_resp / $lower_bound with joint
    # posteriors (see the end of fit_ml), so a conditional fit must have its
    # measurement-only posteriors restored before the new structural model is
    # estimated. One E-step on the frozen parameters recovers them exactly.
    fit$sm <- NULL
    e_res <- e_step(fit, fit$data, NULL)
    fit$log_resp    <- e_res$log_resp
    fit$lower_bound <- e_res$log_prob_norm
  }

  # For an unconditional fit these equal $metrics; for a formerly conditional
  # fit they are recomputed from the just-restored measurement-only posteriors.
  if (is.null(fit$step1_metrics))
    fit$step1_metrics <- .step1_metrics(fit)

  fit
}

# Align user-supplied structural data with the rows the model was actually
# fit on. Cases with no observed indicator are removed before estimation
# (see .empty_rows), so external variables supplied for the original data
# must be subset the same way.
.align_structural_rows <- function(Y, fit, arg_name) {
  md     <- fit$missing_data
  n_kept <- nrow(fit$data)
  n_in   <- md$n_input_rows %||% n_kept
  empty  <- md$empty_rows

  if (nrow(Y) == n_in && length(empty) > 0L) {
    terms_attr <- attr(Y, "covariate_terms")
    Y <- Y[-empty, , drop = FALSE]
    attr(Y, "covariate_terms") <- terms_attr
    message(sprintf(
      paste0("%d case(s) had been removed at fitting time for having no ",
             "observed indicator; the matching rows of `%s` were dropped."),
      length(empty), arg_name))
  } else if (nrow(Y) != n_kept) {
    stop(sprintf(
      paste0("`%s` has %d rows, but the model was fit to %d cases",
             "%s. Supply one row per case of the original data."),
      arg_name, nrow(Y), n_kept,
      if (length(empty) > 0L)
        sprintf(" (%d supplied originally, %d removed for having no observed indicator)",
                n_in, length(empty))
      else ""), call. = FALSE)
  }

  Y
}

# Resolve a one-sided formula or a character vector of column names against a
# `data` argument, returning the named columns as a data frame. This is the only
# thing the `data =` form adds: both verbs rejoin their existing path the moment
# it returns, so nothing downstream knows which form the user typed.
.columns_from_data <- function(spec, data, arg_name) {
  if (!is.data.frame(data)) {
    coerced <- try(as.data.frame(data, stringsAsFactors = FALSE),
                   silent = TRUE)
    if (inherits(coerced, "try-error"))
      stop("`data` must be a data frame, or something coercible to one.",
           call. = FALSE)
    data <- coerced
  }

  if (is.character(spec)) {
    absent <- setdiff(spec, names(data))
    if (length(absent))
      stop(sprintf("`%s` names %s not found in `data`: %s.", arg_name,
                   if (length(absent) > 1L) "columns" else "a column",
                   paste(absent, collapse = ", ")), call. = FALSE)
    return(data[, spec, drop = FALSE])
  }

  if (!inherits(spec, "formula") || length(spec) != 2L)
    stop(sprintf(paste0("When `data` is supplied, `%s` must be a one-sided ",
                        "formula (~ x + y) or a character vector of column ",
                        "names."), arg_name), call. = FALSE)

  absent <- setdiff(all.vars(spec), names(data))
  if (length(absent))
    stop(sprintf("`%s` names %s not found in `data`: %s.", arg_name,
                 if (length(absent) > 1L) "columns" else "a column",
                 paste(absent, collapse = ", ")), call. = FALSE)

  # na.action = na.pass: a case missing a structural variable is kept (it is
  # completed at step 3, see the details of add_covariates()), and dropping
  # rows here would break the one-row-per-case alignment that
  # .align_structural_rows() goes on to check.
  mf <- stats::model.frame(spec, data = data, na.action = stats::na.pass)
  attr(mf, "terms") <- NULL
  mf
}

# Whether a one-sided formula names an interaction (order > 1), the one shape
# `.columns_from_data()`'s model.frame() cannot expand: it returns the
# variables referenced, not the columns a model actually needs, so `a:b` or
# `a*b` silently fits as `a + b` under that path. A plain main-effects formula
# is routed to `.columns_from_data()` unchanged rather than through
# `.covariate_matrix_from_formula()` below, because the two disagree on a
# dummy's column name (`sexo.M` from `prepare_covariates()`'s own convention
# vs. `sexoM` from `model.matrix()`'s), and there is no reason to disturb that
# naming, or the tests and saved output that depend on it, for a formula shape
# that already worked.
.formula_has_interaction <- function(spec) any(attr(stats::terms(spec), "order") > 1L)

# Resolve a one-sided formula against `data` into a numeric design matrix, for
# a *covariate* spec specifically -- never for `outcome` itself, which needs
# the raw variable (see .columns_from_data() above) rather than one already
# dummy-coded, since a categorical outcome's own level structure is what
# .build_outcome_spec() inspects next.
#
# model.matrix() is what .columns_from_data()'s model.frame() cannot give: a
# factor's k-1 dummies and an interaction's several columns collapse to one
# term via `assign`, which is exactly the bookkeeping a hand-built
# model.matrix(~ a * b)[, -1] has no way to carry (see fit_mixture()'s
# `@param group` documentation, which is what tells a user to hand-build one
# in the first place). `na.action = na.pass` at the model.frame stage, kept
# through to model.matrix() via the terms object built from that frame, is
# what lets a case missing a structural covariate stay in the data instead of
# being dropped -- .align_structural_rows() expects one row per case.
.covariate_matrix_from_formula <- function(spec, data, arg_name) {
  absent <- setdiff(all.vars(spec), names(data))
  if (length(absent))
    stop(sprintf("`%s` names %s not found in `data`: %s.", arg_name,
                 if (length(absent) > 1L) "columns" else "a column",
                 paste(absent, collapse = ", ")), call. = FALSE)

  mf   <- stats::model.frame(spec, data = data, na.action = stats::na.pass)
  tt   <- stats::terms(mf)
  mm   <- stats::model.matrix(tt, data = mf)
  asg  <- attr(mm, "assign")
  keep <- asg != 0L                       # drop the intercept column
  out  <- mm[, keep, drop = FALSE]
  attr(out, "covariate_terms") <- attr(tt, "term.labels")[asg[keep]]
  out
}

# Guard the case of a formula with no `data` to resolve it against, which is the
# likeliest way to mistype the new form.
.check_data_form <- function(spec, data, arg_name) {
  if (is.null(data) && inherits(spec, "formula"))
    stop(sprintf("`%s` is a formula, so `data` must be supplied as well.",
                 arg_name), call. = FALSE)
  invisible(NULL)
}

# `steps` on add_covariates()/add_outcome(): 3 is the bias-adjusted three-step,
# 2 the two-step of Bakk and Kuha (2018). `correction` is a property of the
# third step only, so naming one under `steps = 2` is refused rather than
# ignored.
.check_steps <- function(steps, corr_set) {
  if (length(steps) > 1L) steps <- steps[1L]
  if (!steps %in% c(2, 3))
    stop("`steps` must be 3 (bias-adjusted three-step) or 2 (two-step).",
         call. = FALSE)
  if (steps == 2 && corr_set)
    stop("`correction` applies to the three-step only; a two-step fit has no ",
         "classification step to correct. Drop `correction` or set ",
         "`steps = 3`.", call. = FALSE)
  as.integer(steps)
}

# Shared execution: attach the structural model and run steps 2-3 only.
.add_structural <- function(fit, Y_use, engine, correction, se, max_iter,
                            assignment = "proportional",
                            moderated = integer(0), steps = 3L,
                            variances = "equal") {
  fit$sm         <- build_emission(engine, n_components = fit$n_components,
                                   moderated = moderated, variances = variances)
  # The structural model is built here rather than in fit_mixture_internal(), so
  # it has to be handed the fit's prior strengths on the way past; without this
  # `bayes_constants` would apply on the one-call three-step path and silently
  # lapse on this one.
  fit$sm         <- .attach_bayes_constants(
    fit$sm, .resolve_bayes_constants(fit$bayes_constants))
  fit            <- .mirror_design_onto_sm(fit)
  fit$n_steps    <- steps
  fit$correction <- if (steps == 2L) "none" else correction
  # Kept on the fit so a saved object still says which assignment rule produced
  # the correction it reports.
  fit$assignment <- assignment
  # The structural data, as fit_mixture() keeps it on a one-call fit.
  fit$Y          <- Y_use

  fit <- .apply_structural_steps(fit, X = fit$data, Y = Y_use, n_steps = steps,
                                 correction = correction, max_iter = max_iter,
                                 se = se, assignment = assignment)
  # order_by_size = FALSE: re-sorting here could relabel classes relative to
  # the unconditional fit the user has already inspected and reported.
  .finalize_model_state(fit, X = fit$data, Y = Y_use, order_by_size = FALSE)
}

#' Examine Predictors of Class Membership on a Fitted Model
#'
#' @description
#' Takes the latent class model you have already chosen and relates covariates
#' to class membership with the bias-adjusted three-step approach (Vermunt,
#' 2010). The measurement model is reused exactly as fitted — no re-estimation,
#' and no risk of landing on a different solution — so this is both faster and
#' conceptually cleaner than re-specifying the model with `predictors`.
#'
#' @param fit A fitted unconditional model from [fit_mixture()] (or
#'   [fit_rmlca()], [fit_lcga()], [fit_gmm()]).
#' @param predictors Covariates that predict class membership: a data frame,
#'   matrix, or single vector/factor. Factors are dummy-coded with the first
#'   level as reference. Must have one row per case of the data the model was
#'   fit to.
#' @param correction Bias correction for the third step: `"ML"` (default;
#'   Vermunt, 2010), `"BCH"`, or `"none"`. Three-step only; an error with
#'   `steps = 2`.
#' @param steps `3` (default) for the bias-adjusted three-step, or `2` for
#'   the two-step estimator of Bakk and Kuha (2018): `fit`'s measurement
#'   model is held fixed and the class-membership regression is estimated by
#'   maximising the full likelihood, with every case's class probabilities
#'   recomputed under the joint model at each iteration. No classification
#'   step, no correction. `fit` is exactly the step-1 estimate the two-step
#'   starts from. See `n_steps` in [fit_mixture()] for what the two-step is
#'   and is not. Its standard errors carry the step-1 uncertainty (Bakk and
#'   Kuha, 2018, eq. 5); see [covariate_se].
#' @param se Standard-error estimator for the structural coefficients:
#'   `"corrected"` (default), `"robust"`, or `"hessian"`. What each is under
#'   the three-step and under the two-step is in [covariate_se].
#' @param assignment How step 1's posteriors are turned into the assigned-class
#'   variable whose classification error the correction inverts.
#'   `"proportional"` (default) gives every case a weight in every class equal
#'   to its posterior probability; `"modal"` assigns each case to its most
#'   likely class outright. The default follows Bakk, Tekle and Vermunt (2013),
#'   who compared the two rules across 54 simulation conditions and found
#'   proportional at least as accurate everywhere and clearly better when the
#'   classes are poorly separated. Use `"modal"` when reproducing an analysis
#'   whose classes were assigned that way. Three-step only.
#' @param max_iter Maximum iterations for the structural estimation.
#' @param data Optional data frame to take the covariates from, in which case
#'   `predictors` may be a one-sided formula (`~ age + sex`, or `~ age * sex`
#'   for an interaction) or a vector of column names instead of the columns
#'   themselves. A formula's terms -- a factor's dummies, an interaction's
#'   several columns -- are recognised as one term by the omnibus Wald test in
#'   [analytical_wald_test()].
#' @param ... Currently unused.
#'
#' @details
#' A case missing a predictor is retained, not deleted: the missing value is
#' completed under the class-invariant Gaussian marginal of the predictors
#' (Sterba, 2014), so the analysis keeps its full N. An analysis that listwise
#' deletes them is fitted to fewer cases; check the reported N before comparing
#' coefficients with a published set.
#'
#' `se = "corrected"` (the default) is the Bakk, Oberski and Vermunt (2014)
#' estimator, which propagates the uncertainty in the step-1 estimates as well as
#' the step-3 sampling variability. `se = "robust"` reports only the latter; use
#' it when the assigned classes are to be treated as given.
#'
#' **Why a three-step function at all, when the measurement model could just be
#' refit with `predictors` supplied to `fit_mixture()` directly.** Jiang,
#' Elliott, Sammel and Wang (2016) name the failure mode of that one-step
#' alternative: when a covariate participates in forming the classes, the
#' joint models "tend to have a high chance of artificially creating
#' spurious mixture components to enhance predictive accuracy for the sample
#' that is used to derive the model," and the apparent gain does not survive
#' an independent validation sample. Their simulation also found the
#' deviance-based criteria (AIC, BIC) more reliable for choosing the number of
#' classes than predictive ones, which over-selected components in pursuit of
#' fit. `add_covariates()`'s two-stage design -- fit the measurement model
#' first, decide on classes, only then look at what predicts them -- keeps
#' that choice from being made by a covariate that happens to correlate with
#' the outcome.
#'

#' @return A `mixture_model` with the class-membership regression attached.
#'   Use [summary()] for odds ratios and omnibus tests; [coef()],
#'   [confint()], and [wald_omnibus_test()] also apply.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis}, \emph{18}(4), 450–469.
#' \doi{10.1093/pan/mpq025}
#'
#' Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#' assignments to external variables: Standard errors for correct inference.
#' \emph{Political Analysis}, \emph{22}(4), 520–540. \doi{10.1093/pan/mpu003}
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis}, \emph{12}(1), 3–27. \doi{10.1093/pan/mph001}
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#' between latent class membership and external variables using bias-adjusted
#' three-step approaches. \emph{Sociological Methodology}, \emph{43}(1),
#' 272–311. \doi{10.1177/0081175012470644}
#'
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#' classes and external variables. \emph{Psychometrika}, \emph{83}(4),
#' 871–892. \doi{10.1007/s11336-017-9592-7}
#'
#' Jiang, Y., Elliott, M. R., Sammel, M. D., & Wang, N. (2016). Joint modeling
#' of cross-sectional health outcomes and longitudinal predictors via mixtures
#' of latent classes. \emph{Statistics and Its Interface}, \emph{9}, 183–201.
#'
#' @seealso [add_outcome()] for distal outcomes; [fit_mixture()] to fit the
#'   unconditional model.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' age   <- rnorm(100)
#' fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
#' fit_cov <- add_covariates(fit, age)
#' summary(fit_cov)
#'
#' # The same covariate named in a formula against its data frame
#' df <- data.frame(age = age, sex = rbinom(100, 1, 0.5))
#' fit_cov2 <- add_covariates(fit, ~ age + sex, data = df)
#'
#' @export
add_covariates <- function(fit, predictors,
                           correction = c("ML", "BCH", "none"),
                           steps = c(3, 2),
                           se = c("corrected", "robust", "hessian"),
                           assignment = c("proportional", "modal"),
                           max_iter = 1000, data = NULL, ...) {
  corr_set        <- !missing(correction)
  correction      <- match.arg(correction)
  steps           <- .check_steps(steps, corr_set)
  se              <- match.arg(se)
  assignment      <- match.arg(assignment)
  predictors_expr <- substitute(predictors)

  if (missing(predictors) || is.null(predictors))
    stop("`predictors` is required: the covariates that predict class ",
         "membership.", call. = FALSE)

  .check_data_form(predictors, data, "predictors")

  fit <- .check_stepwise_fit(fit, "add_covariates")

  if (!corr_set && steps == 3L)
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))

  if (!is.null(data) && inherits(predictors, "formula") &&
      .formula_has_interaction(predictors)) {
    predictors      <- .covariate_matrix_from_formula(predictors, data, "predictors")
    predictors_expr <- NULL
  } else if (!is.null(data) &&
            (inherits(predictors, "formula") || is.character(predictors))) {
    predictors      <- .columns_from_data(predictors, data, "predictors")
    predictors_expr <- NULL
  }

  Y_use <- prepare_covariates(
    .as_named_covariates(predictors, predictors_expr, "predictor"))
  Y_use <- .align_structural_rows(Y_use, fit, "predictors")

  .add_structural(fit, Y_use, "predict_class", correction, se, max_iter,
                  assignment = assignment, steps = steps)
}

#' Examine a Distal Outcome on a Fitted Model
#'
#' @description
#' Takes the latent class model you have already chosen and relates the classes
#' to a distal outcome with the bias-adjusted three-step approach. The
#' measurement model is reused exactly as fitted — no re-estimation, and no
#' risk of landing on a different solution.
#'
#' @param fit A fitted unconditional model from [fit_mixture()] (or
#'   [fit_rmlca()], [fit_lcga()], [fit_gmm()]).
#' @param outcome The distal outcome: a numeric vector (continuous) or a
#'   factor/character/integer vector (categorical). Must have one value per
#'   case of the data the model was fit to.
#' @param covariates Optional covariates that adjust the outcome.
#' @param outcome_type One of `"auto"` (default; inferred from `outcome`),
#'   `"continuous"`, or `"categorical"`.
#' @param slopes When `covariates` are supplied, whether their effect is
#'   `"pooled"` (default; one slope shared across classes), `"class_specific"`
#'   (every covariate gets its own slope per class), or a character vector of
#'   covariate names (or a one-sided formula naming them, e.g. `~ loc1 +
#'   loc2`) giving a slope per class to just those covariates while the rest
#'   stay pooled. The last form -- letting the class moderate some covariates
#'   while adjusting for others -- is continuous-outcome only.
#' @param predictors Optional covariates that predict class membership, in any
#'   form [add_covariates()] accepts. When supplied, the class-membership
#'   regression and the outcome are estimated together in one step-3 model:
#'   each case's contribution is \eqn{\sum_k P(k \mid z) f(y \mid k, x)} times
#'   the classification-error term, so the class probabilities the outcome is
#'   related to are the covariate-specific ones rather than the overall class
#'   sizes (Vermunt, 2010). `outcome` may then name several distal outcomes (a
#'   data frame, or a formula naming several columns of `data`); each is
#'   specified as it would be on its own and all share `covariates`. Estimated
#'   separately, the outcomes and the predictors would each answer a
#'   different question from different class probabilities; estimated jointly,
#'   the paths from the predictors to the classes, from the classes to each
#'   outcome, and from the covariates to each outcome are adjusted for one
#'   another. Available with `correction = "ML"` (the default here) or
#'   `"BCH"`, and pooled slopes. Under ML the standard errors are those of the
#'   joint step-3 log-likelihood (`se = "hessian"` for the inverse observed
#'   information, otherwise the sandwich); under BCH the model is fitted to
#'   the BCH-weighted log-likelihood and its standard errors are always the
#'   sandwich clustered on the case. Both treat the step-1 estimates as known.
#' @param variances For a continuous outcome with `covariates`: `"equal"`
#'   (default; one residual variance shared by the classes) or
#'   `"class_specific"` (one per class). A continuous outcome without
#'   covariates always has one variance per class.
#' @param correction Bias correction for the third step: `"auto"` (default)
#'   picks `"BCH"` for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"`
#'   for categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.
#'   Three-step only; an error with `steps = 2`.
#' @param steps `3` (default) for the bias-adjusted three-step, or `2` for
#'   the two-step estimator of Bakk and Kuha (2018): `fit`'s measurement
#'   model and class sizes are held fixed and the outcome model is estimated
#'   by maximising the full likelihood, every case's class probabilities
#'   recomputed under the joint model at each iteration. No classification
#'   step, no correction. See `n_steps` in [fit_mixture()]; for a distal
#'   outcome the two-step's standard errors do not yet carry the step-1
#'   uncertainty, and the printed output says so.
#' @param se Standard-error estimator passed on to the third step:
#'   `"corrected"` (default), `"robust"`, or `"hessian"`. It governs the
#'   covariate part of the third step. A continuous distal outcome under
#'   `correction = "BCH"` always reports a sandwich clustered on the case,
#'   whatever this is set to, and so does a BCH model with `predictors`: the
#'   expanded data set carries one weighted
#'   record per class per case, so a case-clustered sandwich is the only
#'   estimator that prices the information the correction gives up.
#' @param assignment How step 1's posteriors are turned into the assigned-class
#'   variable whose classification error the correction inverts.
#'   `"proportional"` (default) gives every case a weight in every class equal
#'   to its posterior probability; `"modal"` assigns each case to its most
#'   likely class outright. The default follows Bakk, Tekle and Vermunt (2013),
#'   who compared the two rules across 54 simulation conditions and found
#'   proportional at least as accurate everywhere and clearly better when the
#'   classes are poorly separated. Use `"modal"` when reproducing an analysis
#'   whose classes were assigned that way.
#' @param max_iter Maximum iterations for the step-3 estimation.
#' @param data Optional data frame to take the variables from, in which case
#'   `outcome` may be a one-sided formula naming one column (`~ bmi`), and
#'   `covariates` a one-sided formula (interactions included, e.g. `~ age *
#'   sex`) or a vector of column names.
#' @param ... Currently unused.
#'
#' @return A `mixture_model` with the distal-outcome model attached. Use
#'   [summary()] for class-specific means or probabilities and their tests, and
#'   [outcome_contrasts()] for which classes differ from which, rather than
#'   whether any of them do.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis}, \emph{18}(4), 450–469.
#' \doi{10.1093/pan/mpq025}
#'
#' Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
#' modeling with continuous distal outcomes. \emph{Structural Equation
#' Modeling}, \emph{23}(1), 20–31. \doi{10.1080/10705511.2014.955104}
#'
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis}, \emph{12}(1), 3–27. \doi{10.1093/pan/mph001}
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#' between latent class membership and external variables using bias-adjusted
#' three-step approaches. \emph{Sociological Methodology}, \emph{43}(1),
#' 272–311. \doi{10.1177/0081175012470644}
#'
#' @seealso [outcome_contrasts()] for class-vs-class differences on the
#'   outcome; [add_covariates()] for predictors of class membership;
#'   [fit_mixture()] to fit the unconditional model.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' bmi   <- rnorm(100, mean = 25)
#' fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
#' fit_out <- add_outcome(fit, bmi)
#' summary(fit_out)
#'
#' # The same outcome named in a formula against its data frame
#' df <- data.frame(bmi = bmi)
#' fit_out2 <- add_outcome(fit, ~ bmi, data = df)
#'
#' \dontrun{
#' # Class moderates level-of-care while age and gender are only adjusted for:
#' # a mix of "class_specific" and "pooled" in one model.
#' fit_mod <- add_outcome(fit, cannabis_days,
#'                        covariates = data.frame(loc1, loc2, loc3, age, gender),
#'                        slopes = c("loc1", "loc2", "loc3"))
#' }
#'
#' @export
add_outcome <- function(fit, outcome, covariates = NULL,
                        outcome_type = c("auto", "continuous", "categorical"),
                        slopes = "pooled", predictors = NULL,
                        variances = c("equal", "class_specific"),
                        correction = c("auto", "BCH", "ML", "none"),
                        steps = c(3, 2),
                        se = c("corrected", "robust", "hessian"),
                        assignment = c("proportional", "modal"),
                        max_iter = 1000, data = NULL, ...) {
  outcome_type <- match.arg(outcome_type)
  corr_set     <- !missing(correction)
  correction   <- match.arg(correction)
  steps        <- .check_steps(steps, corr_set)
  se           <- match.arg(se)
  assignment   <- match.arg(assignment)
  variances    <- match.arg(variances)
  cov_expr     <- substitute(covariates)
  predictors_expr <- substitute(predictors)

  if (missing(outcome) || is.null(outcome))
    stop("`outcome` is required: the distal outcome to relate to the classes.",
         call. = FALSE)

  .check_data_form(outcome, data, "outcome")
  .check_data_form(covariates, data, "covariates")

  fit <- .check_stepwise_fit(fit, "add_outcome")

  # A character `outcome` is a categorical outcome, not a column name — hence
  # the formula-only test here, where `predictors` and `covariates` take either.
  if (!is.null(data) && inherits(outcome, "formula")) {
    outcome <- .columns_from_data(outcome, data, "outcome")
    if (ncol(outcome) != 1L && is.null(predictors))
      stop(sprintf(paste0("`outcome` must name exactly one distal outcome, ",
                          "but %d were named (%s). Fit them one at a time."),
                   ncol(outcome), paste(names(outcome), collapse = ", ")),
           call. = FALSE)
  }
  if (!is.null(data) && inherits(covariates, "formula") &&
      .formula_has_interaction(covariates)) {
    covariates <- .covariate_matrix_from_formula(covariates, data, "covariates")
    cov_expr   <- NULL
  } else if (!is.null(data) &&
            (inherits(covariates, "formula") || is.character(covariates))) {
    covariates <- .columns_from_data(covariates, data, "covariates")
    cov_expr   <- NULL
  }

  if (!is.null(predictors))
    return(.add_joint_outcome(fit, outcome, covariates, cov_expr, predictors,
                              predictors_expr, data, outcome_type, slopes,
                              variances, correction, corr_set, steps, se,
                              assignment, max_iter))

  spec <- .build_outcome_spec(outcome, covariates, outcome_type, slopes,
                              cov_expr)
  spec <- .apply_outcome_variances(spec, variances)

  if (steps == 2L) {
    correction <- "none"
  } else if (correction == "auto") {
    correction <- if (startsWith(spec$engine, "categorical")) "ML" else "BCH"
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))
  }

  Y_use <- .align_structural_rows(spec$Y, fit, "outcome")

  .add_structural(fit, Y_use, spec$engine, correction, se, max_iter,
                  assignment = assignment, moderated = spec$moderated,
                  steps = steps, variances = spec$variances)
}

# `variances = "class_specific"` exists only where the residual variance is a
# choice: a continuous outcome adjusted for covariates. A continuous outcome
# without covariates already has one variance per class, and a categorical one
# has none.
.apply_outcome_variances <- function(spec, variances) {
  if (identical(variances, "class_specific") &&
      spec$engine != "continuous_outcome_adjusted")
    stop('`variances = "class_specific"` applies to a continuous outcome with ',
         '`covariates` and pooled or partly pooled `slopes`. A continuous ',
         'outcome without covariates already has one variance per class.',
         call. = FALSE)
  spec$variances <- variances
  spec
}

# add_outcome() with `predictors`: the class-membership regression and every
# distal outcome in one step-3 model. The outcomes are not independent pieces
# here the way they are without predictors -- each one's likelihood carries
# the class probabilities P(k | z) that the predictors now model -- so they are
# estimated together, in one nested structural model: a covariate block
# followed by one block per outcome, all fitted by the same ML-corrected EM.
.add_joint_outcome <- function(fit, outcome, covariates, cov_expr, predictors,
                               predictors_expr, data, outcome_type, slopes,
                               variances, correction, corr_set, steps, se,
                               assignment, max_iter) {
  if (steps == 2L)
    stop("`predictors` with a distal outcome is available for the three-step ",
         "(`steps = 3`) only.", call. = FALSE)
  if (correction == "none")
    stop(paste0(
      "`correction = \"none\"` is not yet available when `predictors` and an ",
      "outcome are estimated together; use `correction = \"ML\"` or \"BCH\"."),
      call. = FALSE)
  if (correction == "auto") {
    correction <- "ML"
    message("Using 'ML' bias correction (set `correction` to override).")
  }

  .check_data_form(predictors, data, "predictors")
  if (!is.null(data) && inherits(predictors, "formula") &&
      .formula_has_interaction(predictors)) {
    predictors      <- .covariate_matrix_from_formula(predictors, data, "predictors")
    predictors_expr <- NULL
  } else if (!is.null(data) &&
             (inherits(predictors, "formula") || is.character(predictors))) {
    predictors      <- .columns_from_data(predictors, data, "predictors")
    predictors_expr <- NULL
  }
  Z <- prepare_covariates(
    .as_named_covariates(predictors, predictors_expr, "predictor"))

  # Several outcomes may be named here. Each is specified exactly as a single
  # add_outcome() call would specify it, and all share `covariates`.
  outs <- if (is.data.frame(outcome)) as.list(outcome)
          else if (is.matrix(outcome))
            stats::setNames(lapply(seq_len(ncol(outcome)), function(j) outcome[, j]),
                            colnames(outcome) %||% paste0("outcome", seq_len(ncol(outcome))))
          else list(outcome)
  specs <- lapply(seq_along(outs), function(j) {
    o <- outs[[j]]
    if (length(outs) > 1L) {
      o <- data.frame(o)
      names(o) <- names(outs)[j]
    }
    s <- .build_outcome_spec(o, covariates, outcome_type, slopes, cov_expr)
    if (!s$engine %in% c("continuous_outcome", "continuous_outcome_adjusted",
                         "categorical_outcome", "categorical_outcome_adjusted"))
      stop("With `predictors`, the outcome's covariate slopes must be pooled ",
           '(`slopes = "pooled"`, or names for a continuous outcome); ',
           '`slopes = "class_specific"` is not available in a joint model.',
           call. = FALSE)
    # A categorical outcome has no residual variance, so `variances` speaks
    # only to the continuous outcomes among several.
    if (startsWith(s$engine, "categorical") && length(outs) > 1L) {
      s$variances <- "equal"
      s
    } else .apply_outcome_variances(s, variances)
  })

  blocks <- list(predictor = list(model = "predict_class", n_columns = ncol(Z)))
  nm     <- c("distal", if (length(specs) > 1L)
    paste0("distal", seq_len(length(specs) - 1L) + 1L))
  for (j in seq_along(specs))
    blocks[[nm[j]]] <- list(model = specs[[j]]$engine,
                            n_columns = ncol(specs[[j]]$Y),
                            moderated = specs[[j]]$moderated,
                            variances = specs[[j]]$variances)
  Y_use <- do.call(.cbind_covariates, c(list(Z), lapply(specs, `[[`, "Y")))
  Y_use <- .align_structural_rows(Y_use, fit, "outcome")

  .add_structural(fit, Y_use, blocks, correction, se, max_iter,
                  assignment = assignment, steps = steps)
}
