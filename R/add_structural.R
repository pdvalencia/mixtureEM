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

# Shared execution: attach the structural model and run steps 2-3 only.
.add_structural <- function(fit, Y_use, engine, correction, se, max_iter) {
  fit$sm         <- build_emission(engine, n_components = fit$n_components)
  fit            <- .mirror_design_onto_sm(fit)
  fit$n_steps    <- 3L
  fit$correction <- correction

  fit <- .apply_structural_steps(fit, X = fit$data, Y = Y_use, n_steps = 3L,
                                 correction = correction, max_iter = max_iter,
                                 se = se)
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
#' @param max_iter Maximum iterations for the step-3 estimation.
#' @param ... Currently unused.
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
#' @seealso [add_outcome()] for distal outcomes; [fit_mixture()] to fit the
#'   unconditional model.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' age   <- rnorm(100)
#' fit   <- fit_mixture(items, n_classes = 2)
#' fit_cov <- add_covariates(fit, age)
#' summary(fit_cov)
#'
#' @export
add_covariates <- function(fit, predictors,
                           correction = c("ML", "BCH", "none"),
                           se = c("corrected", "robust", "hessian"),
                           max_iter = 1000, ...) {
  corr_set        <- !missing(correction)
  correction      <- match.arg(correction)
  se              <- match.arg(se)
  predictors_expr <- substitute(predictors)

  if (missing(predictors) || is.null(predictors))
    stop("`predictors` is required: the covariates that predict class ",
         "membership.", call. = FALSE)

  fit <- .check_stepwise_fit(fit, "add_covariates")

  if (!corr_set)
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))

  Y_use <- prepare_covariates(
    .as_named_covariates(predictors, predictors_expr, "predictor"))
  Y_use <- .align_structural_rows(Y_use, fit, "predictors")

  .add_structural(fit, Y_use, "predict_class", correction, se, max_iter)
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
#'   `"pooled"` (one slope shared across classes) or `"class_specific"`.
#' @param correction Bias correction for the third step: `"auto"` (default)
#'   picks `"BCH"` for continuous outcomes (Bakk & Vermunt, 2016) and `"ML"`
#'   for categorical outcomes; or set `"BCH"`, `"ML"`, `"none"` directly.
#' @param se Standard-error estimator passed on to the third step:
#'   `"corrected"` (default), `"robust"`, or `"hessian"`.
#' @param max_iter Maximum iterations for the step-3 estimation.
#' @param ... Currently unused.
#'
#' @return A `mixture_model` with the distal-outcome model attached. Use
#'   [summary()] for class-specific means or probabilities and their tests.
#'
#' @references
#' Bakk, Z., & Vermunt, J. K. (2016). Robustness of stepwise latent class
#' modeling with continuous distal outcomes. \emph{Structural Equation
#' Modeling}, \emph{23}(1), 20–31. \doi{10.1080/10705511.2014.955104}
#'
#' @seealso [add_covariates()] for predictors of class membership;
#'   [fit_mixture()] to fit the unconditional model.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' bmi   <- rnorm(100, mean = 25)
#' fit   <- fit_mixture(items, n_classes = 2)
#' fit_out <- add_outcome(fit, bmi)
#' summary(fit_out)
#'
#' @export
add_outcome <- function(fit, outcome, covariates = NULL,
                        outcome_type = c("auto", "continuous", "categorical"),
                        slopes = c("pooled", "class_specific"),
                        correction = c("auto", "BCH", "ML", "none"),
                        se = c("corrected", "robust", "hessian"),
                        max_iter = 1000, ...) {
  outcome_type <- match.arg(outcome_type)
  slopes       <- match.arg(slopes)
  correction   <- match.arg(correction)
  se           <- match.arg(se)
  cov_expr     <- substitute(covariates)

  if (missing(outcome) || is.null(outcome))
    stop("`outcome` is required: the distal outcome to relate to the classes.",
         call. = FALSE)

  fit <- .check_stepwise_fit(fit, "add_outcome")

  spec <- .build_outcome_spec(outcome, covariates, outcome_type, slopes,
                              cov_expr)

  if (correction == "auto") {
    correction <- if (startsWith(spec$engine, "categorical")) "ML" else "BCH"
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))
  }

  Y_use <- .align_structural_rows(spec$Y, fit, "outcome")

  .add_structural(fit, Y_use, spec$engine, correction, se, max_iter)
}
