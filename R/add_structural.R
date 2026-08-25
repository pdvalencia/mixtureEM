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

# Shared execution: attach the structural model and run steps 2-3 only.
.add_structural <- function(fit, Y_use, engine, correction, se, max_iter,
                            assignment = "proportional",
                            moderated = integer(0)) {
  fit$sm         <- build_emission(engine, n_components = fit$n_components,
                                   moderated = moderated)
  # The structural model is built here rather than in fit_mixture_internal(), so
  # it has to be handed the fit's prior strengths on the way past; without this
  # `bayes_constants` would apply on the one-call three-step path and silently
  # lapse on this one.
  fit$sm         <- .attach_bayes_constants(
    fit$sm, .resolve_bayes_constants(fit$bayes_constants))
  fit            <- .mirror_design_onto_sm(fit)
  fit$n_steps    <- 3L
  fit$correction <- correction
  # Kept on the fit so a saved object still says which assignment rule produced
  # the correction it reports.
  fit$assignment <- assignment

  fit <- .apply_structural_steps(fit, X = fit$data, Y = Y_use, n_steps = 3L,
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
#'   Vermunt, 2010), `"BCH"`, or `"none"`.
#' @param se Standard-error estimator passed on to the third step:
#'   `"corrected"` (default), `"robust"`, or `"hessian"`.
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
#' it when reproducing an analysis whose standard errors were computed that way.
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
                           se = c("corrected", "robust", "hessian"),
                           assignment = c("proportional", "modal"),
                           max_iter = 1000, data = NULL, ...) {
  corr_set        <- !missing(correction)
  correction      <- match.arg(correction)
  se              <- match.arg(se)
  assignment      <- match.arg(assignment)
  predictors_expr <- substitute(predictors)

  if (missing(predictors) || is.null(predictors))
    stop("`predictors` is required: the covariates that predict class ",
         "membership.", call. = FALSE)

  .check_data_form(predictors, data, "predictors")

  fit <- .check_stepwise_fit(fit, "add_covariates")

  if (!corr_set)
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
                  assignment = assignment)
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
#' @param correction Bias correction for the third step: `"auto"` (default)
#'   picks `"BCH"` for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"`
#'   for categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.
#' @param se Standard-error estimator passed on to the third step:
#'   `"corrected"` (default), `"robust"`, or `"hessian"`. It governs the
#'   covariate part of the third step. A continuous distal outcome under
#'   `correction = "BCH"` always reports a sandwich clustered on the case,
#'   whatever this is set to: the expanded data set carries one weighted
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
                        slopes = "pooled",
                        correction = c("auto", "BCH", "ML", "none"),
                        se = c("corrected", "robust", "hessian"),
                        assignment = c("proportional", "modal"),
                        max_iter = 1000, data = NULL, ...) {
  outcome_type <- match.arg(outcome_type)
  correction   <- match.arg(correction)
  se           <- match.arg(se)
  assignment   <- match.arg(assignment)
  cov_expr     <- substitute(covariates)

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
    if (ncol(outcome) != 1L)
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

  spec <- .build_outcome_spec(outcome, covariates, outcome_type, slopes,
                              cov_expr)

  if (correction == "auto") {
    correction <- if (startsWith(spec$engine, "categorical")) "ML" else "BCH"
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))
  }

  Y_use <- .align_structural_rows(spec$Y, fit, "outcome")

  .add_structural(fit, Y_use, spec$engine, correction, se, max_iter,
                  assignment = assignment, moderated = spec$moderated)
}
