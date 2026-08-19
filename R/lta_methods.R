# ==============================================================================
# Latent transition analysis - extractors, tests, and display methods
# ==============================================================================

# ------------------------------------------------------------------------------
# Extractors
# ------------------------------------------------------------------------------

#' Transition Probability Matrices
#'
#' @description
#' Returns the estimated transition probabilities \eqn{\tau_{s_t | s_{t-1}}}: the
#' probability of occupying each latent status at one occasion given the status
#' held at the previous one. Rows are the origin status and sum to one; the
#' diagonal is the probability of staying put.
#'
#' There is one matrix per pair of adjacent occasions unless the model was fitted
#' with `tau_homogeneous = TRUE`, in which case a single matrix is shared.
#'
#' In a mixture latent Markov model each latent class has its own transitions,
#' so the result gains an outer level indexed by class; `class` picks one out.
#'
#' @param object An object returned by [`fit_lta()`].
#' @param occasion Optional index of a single transition (1 for the move from
#'   occasion 1 to occasion 2, and so on). Omit to get them all.
#' @param class Optional latent class, for a model fitted with `n_classes` > 1
#'   or `mover_stayer = TRUE`. Omit to get every class.
#' @return A matrix, or a named list of matrices when `occasion` is omitted and
#'   the transitions are time-heterogeneous, nested inside a list over classes
#'   when the model has more than one.
#' @export
transition_matrix <- function(object, occasion = NULL, class = NULL) {
  if (!inherits(object, "lta_model"))
    stop("`object` must be a fitted latent transition model.", call. = FALSE)
  labs <- object$longitudinal$time_labels
  K    <- object$n_statuses
  C    <- object$n_classes %||% 1L
  st   <- paste0("Status ", seq_len(K))

  one_class <- function(taus) {
    out <- lapply(taus, function(m) {
      dimnames(m) <- list(from = st, to = st)
      m
    })
    names(out) <- sprintf("%s -> %s", labs[-length(labs)], labs[-1])
    # A shared table is returned as a bare matrix rather than a one-element
    # list, and without an extra attribute: print() already reports that the
    # transitions are held equal, and a stray attribute under the table is
    # noise.
    if (isTRUE(object$tau_homogeneous) && is.null(occasion)) return(out[[1]])
    if (!is.null(occasion)) return(out[[occasion]])
    out
  }

  if (C == 1L) return(one_class(object$tau))
  if (!is.null(class)) return(one_class(object$tau[[class]]))
  out <- lapply(object$tau, one_class)
  names(out) <- .lta_class_labels(object)
  out
}

#' @rdname class_assignments
#'
#' @details
#' For a latent transition model the assignment is of latent *status*, and a
#' status assignment is made at every occasion rather than once per case — the
#' same convention [`status_prevalences()`] and the entropy in
#' [`fit_lta()`]'s `metrics` already follow. Supply `occasion` to work with one
#' occasion at the shape the mixture methods return; omit it for all of them at
#' once. To assign the latent *class* of a mixture latent Markov model, use
#' `object$class_posterior`.
#'
#' @param occasion For an `lta_model`, the index of a single occasion. Omit for
#'   every occasion, which `type = "both"` does not support.
#' @export
class_assignments.lta_model <- function(object,
                                        type = c("modal", "posterior", "both"),
                                        occasion = NULL, ...) {
  type <- match.arg(type)
  gam  <- object$gamma
  labs <- object$longitudinal$time_labels
  st   <- paste0("Status ", seq_len(object$n_statuses))

  if (!is.null(occasion)) {
    if (length(occasion) != 1L || !occasion %in% seq_along(gam))
      stop("`occasion` must be a single index between 1 and ", length(gam),
           ".", call. = FALSE)
    g <- gam[[occasion]]
    modal <- max.col(g, ties.method = "first")
    if (type == "modal") return(modal)
    colnames(g) <- st
    if (type == "posterior") return(g)
    return(data.frame(status      = modal,
                      probability = g[cbind(seq_len(nrow(g)), modal)],
                      g,
                      check.names = FALSE))
  }

  if (type == "both")
    stop("`type = \"both\"` needs a single `occasion`: an LTA case has one ",
         "status per occasion, not one overall.", call. = FALSE)

  if (type == "posterior") {
    out <- lapply(gam, function(g) { colnames(g) <- st; g })
    names(out) <- labs
    return(out)
  }

  out <- vapply(gam, max.col, ties.method = "first",
                FUN.VALUE = integer(nrow(gam[[1]])))
  colnames(out) <- labs
  out
}

# Class labels, naming the stayer where there is one. The stayer is the last
# class by construction; see .lta_tau_allowed().
.lta_class_labels <- function(object) {
  C <- object$n_classes %||% 1L
  labs <- paste("Class", seq_len(C))
  if (isTRUE(object$mover_stayer)) {
    labs[C] <- paste0(labs[C], " (stayer)")
    labs[-C] <- paste0(labs[-C], " (mover)")
  }
  labs
}

#' Latent Status Prevalences by Occasion
#'
#' @description
#' The proportion in each latent status at each occasion. The model-implied
#' prevalences propagate \eqn{\delta} through the transition matrices; the
#' empirical ones average the posterior status probabilities.
#'
#' With more than one latent class the whole-sample prevalences are returned by
#' default - the mixture's own marginal, which is what should be compared with
#' the observed proportions - and `class` picks out a single class's chain.
#'
#' @param object An object returned by [`fit_lta()`].
#' @param type `"model"` (default) or `"posterior"`.
#' @param class Optional latent class, for a model fitted with `n_classes` > 1
#'   or `mover_stayer = TRUE`.
#' @return An occasions-by-statuses matrix.
#' @export
status_prevalences <- function(object, type = c("model", "posterior"),
                               class = NULL) {
  if (!inherits(object, "lta_model"))
    stop("`object` must be a fitted latent transition model.", call. = FALSE)
  type <- match.arg(type)
  C <- object$n_classes %||% 1L
  if (!is.null(class) && C == 1L)
    stop("`class` applies only to a model with more than one latent class.",
         call. = FALSE)

  if (type == "model") {
    if (is.null(class)) return(.lta_bare_prevalences(object$prevalences))
    return(attr(object$prevalences, "by_class")[[class]])
  }

  w <- object$weights_vec
  # A class's empirical prevalence weights each case by how likely it is to
  # belong to that class, which is what makes the per-class curves add back up
  # to the whole-sample one.
  gam <- if (is.null(class)) object$gamma else object$gamma_by_class[[class]]
  wc  <- if (is.null(class)) w else w * object$class_posterior[, class]
  P <- t(vapply(gam, function(g) colSums(g * wc) / sum(wc),
                numeric(object$n_statuses)))
  dimnames(P) <- dimnames(.lta_bare_prevalences(object$prevalences))
  P
}

# ------------------------------------------------------------------------------
# Absolute fit
# ------------------------------------------------------------------------------

#' Likelihood-Ratio Goodness-of-Fit Statistic
#'
#' @description
#' Computes the likelihood-ratio chi-square \eqn{G^2} (also written \eqn{L^2})
#' comparing the observed response-pattern frequencies with those the model
#' implies, together with its degrees of freedom
#' \eqn{df = W - P - 1}, where \eqn{W} is the number of cells in the contingency
#' table formed by crossing every item at every occasion and \eqn{P} the number
#' of free parameters (Collins & Lanza, 2010, sec. 4.3.2 and 7.6).
#'
#' The statistic is defined only for fully categorical indicators observed
#' without missingness. Even then it should be read with care: the table has
#' \eqn{W} cells and is usually extremely sparse, so the chi-square reference
#' distribution is unreliable and the value is best used to compare models rather
#' than to test one in isolation.
#'
#' This is the longitudinal-facing name for [`absolute_fit()`], which computes
#' the same \eqn{G^2} for any categorical mixture model and reports the Pearson
#' \eqn{X^2} and Cressie-Read statistics alongside it. The two are
#' interchangeable; this one returns a plain list for backward compatibility.
#'
#' @param object A model fitted by [`fit_lta()`] or [`fit_rmlca()`].
#' @return A list with `g2`, `df`, `p_value`, `n_cells` and `n_patterns`, or
#'   `NULL` (with a message) when the statistic does not apply.
#' @seealso [`absolute_fit()`], [`bivariate_residuals()`].
#' @export
lta_g2 <- function(object) {
  fit <- absolute_fit(object)
  if (is.null(fit)) return(NULL)
  fit[c("g2", "df", "p_value", "n_cells", "n_patterns")]
}

# Number of response categories for each column of the indicator matrix, or
# NULL if any indicator is not categorical. Derived from the same accessor the
# bivariate residuals use (R/fit_diagnostics.R), so a measurement model that can
# produce a two-way table can produce a full one and vice versa.
.longitudinal_col_levels <- function(mm, n_cols) {
  items <- .categorical_item_probs(mm)
  if (is.null(items) || length(items) != n_cols) return(NULL)
  vapply(items, function(it) length(it$categories), integer(1))
}

# Uniform accessor over the model classes that can have a contingency table.
# `conditional` marks a fit whose case-level probabilities depend on covariates,
# for which no single model-implied table exists.
.nested_fit_info <- function(object) {
  if (inherits(object, "lta_model"))
    return(list(ll = object$loglik, n_params = object$n_params,
                data = object$data, mm = object$mm,
                weights = object$weights_vec, ll_case = object$ll_case,
                conditional = !is.null(object$Z_delta) ||
                              !is.null(object$Z_tau),
                label = "LTA"))
  if (inherits(object, "mixture_model"))
    return(list(ll = object$metrics$ll, n_params = object$metrics$n_params,
                data = object$data, mm = object$mm,
                weights = object$sample_weights,
                ll_case = object$lower_bound,
                conditional = has_covariate(object$sm),
                label = if (inherits(object, "rmlca")) "RMLCA" else "Mixture"))
  stop("Unsupported model object.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# Choosing the number of classes / statuses
# ------------------------------------------------------------------------------

#' Compare Longitudinal Mixture Models Across a Range of Class Counts
#'
#' @description
#' Fits a series of models with increasing numbers of latent classes (RMLCA) or
#' latent statuses (LTA) and tabulates the usual selection criteria, mirroring
#' [`compare_mixtures()`] for the cross-sectional case. Lower AIC, BIC and SABIC
#' are better; entropy summarises how cleanly cases are classified.
#'
#' Selecting the number of latent statuses for an LTA should use the data from
#' every occasion at once rather than a separate cross-sectional analysis per
#' occasion: pooling the repeated measures gives the model more information, so
#' a solution can be supported longitudinally that no single occasion would
#' support on its own (Collins & Lanza, 2010, sec. 7.3.3).
#'
#' Note that likelihood-based tests of \eqn{K} against \eqn{K-1} classes, such
#' as the bootstrap likelihood-ratio test, do not have their usual reference
#' distribution here, which is why only information criteria are reported.
#'
#' For the growth models the comparison should include the one-class solution,
#' which is the ordinary latent growth curve model: it is the benchmark the
#' class solutions have to beat, and reporting it is asked for by name in the
#' GRoLTS reporting checklist (van de Schoot et al., 2017, item 11). It is
#' therefore in the default `k_range` for `"gmm"` and `"lcga"` and not for the
#' other two, where a one-status LTA is not a model anyone reports.
#'
#' The table informs the decision; it does not make it. As Ram and Grimm (2009,
#' p. 571) put it, "there is not a deterministic set of rules to follow when
#' selecting the best model. Rather, model selection is an art — informed by
#' theory, past findings, past experience, and a variety of statistical fit
#' indices."
#'
#' **Reading the `Entropy` column.** Relative entropy describes how cleanly a
#' solution separates its classes. The usual anchors are 0.40, 0.60 and 0.80 for
#' low, medium and high separation (Clark & Muthen, 2009, as reported by Lee et
#' al., 2023, p. 653), and Ram and Grimm (p. 571) suggest preferring the
#' higher-entropy model when choosing among models with similar BIC. Those
#' anchors are on the same normalisation this package uses. Entropy is not
#' evidence for how many classes there are, though, and there is no threshold it
#' has to clear — "there are no set cut-off criteria for deciding whether the
#' entropy is reasonably high" (Jung & Wickrama, 2008, p. 312) — so the package
#' applies none.
#'
#' **Reading the `Unreplicated` column.** `TRUE` means that K's reported maximum
#' was found by exactly one random start; refit those with more starts
#' (`n_init = 100` is the usual next step) before reporting. The per-model
#' warning is suppressed inside this loop, since
#' it would otherwise fire once per K.
#'
#' @param indicators The repeated indicators, in any format accepted by
#'   [`fit_rmlca()`]; for `model = "gmm"` or `"lcga"`, the single repeated
#'   outcome, in any format accepted by [`fit_gmm()`].
#' @param k_range Integer vector of class or status counts to fit. Defaults to
#'   `1:4` for the growth models and `2:4` for the others.
#' @param model `"lta"` (default), `"rmlca"`, `"gmm"` or `"lcga"`.
#' @param times Number of occasions; required for wide input.
#' @param verbose Print progress.
#' @param vlmr Whether to add the Vuong-Lo-Mendell-Rubin test of K against K+1:
#'   `"none"` (the default), `"standard"`, `"robust"` or `"both"`. See
#'   [`compare_mixtures()`] for what the two forms are and why the test is off
#'   by default. It is unavailable for `model = "lta"`, whose latent variable is
#'   a status per occasion rather than one class per case; those rows come back
#'   `NA` with a message.
#' @param ... Further arguments passed to the fitting function, such as
#'   `measurement`, `time_invariance` or `tau_homogeneous` for the categorical
#'   models, `degree`, `time_scores`, `random_effects`, `psi` or `residual` for
#'   the growth models, and `n_init` for any of them. Passing the same
#'   specification to every K is the point: the criteria in the table are only
#'   comparable across models that differ in nothing else. Ram and Grimm (2009,
#'   p. 571) state the same rule for the tests: they "compare models that differ
#'   only in the number of classes ... but are not appropriate for comparing
#'   models that allow for different types of between-class differences".
#'
#' @return An object of class `mixture_comparison`, a list with `fit_table`
#'   (columns `Classes`, `LL`, `Params`, `AIC`, `BIC`, `SABIC`, `Entropy` and
#'   `Unreplicated`), the fitted `models` (named `"K2"`, `"K3"`, ...) and
#'   `best_k`, the class count with the lowest BIC. It indexes exactly as a
#'   plain list; [`plot()`][plot.mixture_comparison] draws the criteria
#'   against K.
#' @references
#' van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., & Vermunt,
#' J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on latent
#' trajectory studies. \emph{Structural Equation Modeling}, \emph{24}(3),
#' 451-467.
#'
#' Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class growth
#' analysis and growth mixture modeling. \emph{Social and Personality Psychology
#' Compass}, \emph{2}(1), 302-317. \doi{10.1111/j.1751-9004.2007.00054.x}
#'
#' Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction to
#' growth mixture models (GMM). In \emph{International Encyclopedia of
#' Education} (4th ed., Vol. 14, pp. 646-655). Elsevier.
#' \doi{10.1016/B978-0-12-818630-5.10076-4}
#'
#' Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
#' identifying differences in longitudinal change among unobserved groups.
#' \emph{International Journal of Behavioral Development}, \emph{33}(6),
#' 565-576. \doi{10.1177/0165025409343765}
#' @export
compare_longitudinal <- function(indicators, k_range = NULL,
                                 model = c("lta", "rmlca", "gmm", "lcga"),
                                 times = NULL, verbose = TRUE,
                                 vlmr = c("none", "standard", "robust",
                                          "both"), ...) {
  model <- match.arg(model)
  vlmr  <- match.arg(vlmr)
  growth <- model %in% c("gmm", "lcga")
  if (is.null(k_range)) k_range <- if (growth) 1:4 else 2:4
  if (any(k_range < 1L))
    stop("`k_range` must be positive.", call. = FALSE)

  if (verbose)
    message(sprintf("Comparing %s models across K = %s ...",
                    toupper(model), paste(range(k_range), collapse = " to ")))

  models <- list()
  rows   <- list()
  for (k in k_range) {
    if (verbose) message(sprintf("  Fitting %d-%s model...", k,
                                 if (model == "lta") "status" else "class"))
    # Every fit here warns on its own about an unreplicated maximum, so a
    # four-K comparison would raise the same warning four times before the
    # table it belongs next to has even been printed. The column below carries
    # the same information, per K, in the place the user is reading. The
    # transition prior's warning is muffled for the same reason, and it is the
    # more insistent of the two: adding statuses divides the sample among more
    # origin rows, so a sweep raises it on every K above the one it starts at.
    # print() on the chosen model still shows it.
    fit <- withCallingHandlers(
      switch(model,
        lta   = fit_lta(indicators, n_statuses = k, times = times, ...),
        rmlca = fit_rmlca(indicators, n_classes = k, times = times, ...),
        gmm   = fit_gmm(indicators, n_classes = k, times = times, ...),
        lcga  = fit_lcga(indicators, n_classes = k, times = times, ...)),
      mixtureEM_replication = function(w) invokeRestart("muffleWarning"),
      mixtureEM_smoothing   = function(w) invokeRestart("muffleWarning"))

    m <- fit$metrics
    models[[paste0("K", k)]] <- fit
    # A one-class model has no classification to be entropic about, so the cell
    # is NA rather than a 1 that would read as perfect separation.
    rows[[length(rows) + 1L]] <- data.frame(
      Classes = k, LL = m$ll, Params = m$n_params,
      AIC = m$aic, BIC = m$bic, SABIC = m$sabic,
      Entropy = if (k == 1L) NA_real_ else m$entropy %||% NA_real_,
      Unreplicated = .is_unreplicated(m))
  }

  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$Classes), , drop = FALSE]
  rownames(tab) <- NULL
  # `lta` is the one model here that does not ride the mixture engine, so
  # .vlmr_pair() declines it by class rather than by a special case.
  if (vlmr != "none") tab <- .vlmr_augment(tab, models, vlmr)
  best <- tab$Classes[which.min(tab$BIC)]
  if (verbose) {
    message("\n=== Model Selection Summary ===")
    print(format(tab, digits = 4))
    flagged <- tab$Classes[which(tab$Unreplicated)]
    if (length(flagged)) {
      # Every K in one call is fitted at the same `n_init`, so the advice is
      # read off the first flagged model rather than repeated per row.
      m1 <- models[[paste0("K", flagged[1])]]$metrics
      message(sprintf("Unreplicated maximum at K = %s - %s.",
                      paste(flagged, collapse = ", "),
                      .replication_advice(m1$n_requested %||% m1$n_starts)))
    }
    message(sprintf("\n-> Best model according to BIC: %d", best))
  }
  # Same class as compare_mixtures() returns, so one plot method serves both.
  out <- list(fit_table = tab, models = models, best_k = best)
  if (vlmr != "none") {
    out$vlmr <- attr(tab, "vlmr_detail")
    attr(out$fit_table, "vlmr_detail") <- NULL
  }
  class(out) <- "mixture_comparison"
  out
}

# ------------------------------------------------------------------------------
# Nested-model tests
# ------------------------------------------------------------------------------

#' Likelihood-Ratio Test for Two Nested Models
#'
#' @description
#' Compares two nested models by the likelihood-ratio difference test,
#' \eqn{-2(\ell_0 - \ell_1)} on \eqn{P_1 - P_0} degrees of freedom. It accepts
#' any pair of fits from [`fit_mixture()`], [`fit_rmlca()`] or [`fit_lta()`],
#' and answers questions of the form "does freeing these parameters buy a
#' significantly better fit?"
#'
#' \itemize{
#'   \item **Measurement invariance across groups** (Collins & Lanza, 2010, sec.
#'     5.8): fit [`fit_mixture()`] with `group_effects = "prevalence"` and with
#'     `"both"`, and compare. This is a cross-sectional test.
#'   \item **Equal prevalences across groups** (sec. 5.11): compare
#'     `group_effects = "none"` against `"prevalence"`.
#'   \item **Measurement invariance across time** (sec. 7.11): fit [`fit_lta()`]
#'     with `measurement_invariance = "full"` and `"none"` and compare.
#'   \item **A time-homogeneous transition matrix** (sec. 7.14): fit
#'     [`fit_lta()`] with `transition_invariance = "full"` and `"none"`.
#' }
#'
#' The models must be nested and fitted to the same data. That is not checked
#' beyond the parameter counts and sample size, so it remains the analyst's
#' responsibility.
#'
#' Because `full` strictly nests `restricted`, its log-likelihood can never
#' be genuinely lower — if it comes out that way here, the `full` model's
#' random-restart search landed on a worse local optimum than the
#' `restricted` model's did, not a real result. A warning is issued in that
#' case; refitting `full` with a larger `n_init` is the usual fix.
#'
#' @param restricted The more constrained model (fewer parameters).
#' @param full The less constrained model.
#' @return A list of class `"lr_test"`.
#' @references
#' Collins, L. M., & Lanza, S. T. (2010). \emph{Latent Class and Latent
#' Transition Analysis: With Applications in the Social, Behavioral, and Health
#' Sciences}. Wiley.
#' @export
lr_test <- function(restricted, full) {
  a <- .nested_fit_info(restricted)
  b <- .nested_fit_info(full)

  if (a$n_params > b$n_params)
    stop("`restricted` has more parameters than `full`; the arguments look ",
         "reversed.", call. = FALSE)
  if (length(a$weights) != length(b$weights))
    stop("The two models were fitted to different numbers of cases.",
         call. = FALSE)

  # A collapsed class variance in either fit invalidates the test in both
  # directions, so it is checked before the sign of the statistic is even looked
  # at. A degenerate fit's log-likelihood is not on the same scale as an
  # admissible one -- it is a spurious optimum sitting on a spike in the
  # likelihood -- so neither a large statistic nor a small one means anything.
  # See R/gaussian_boundary.R.
  # A growth model records the same finding in its own field, since its remedy
  # is structural rather than the variance prior; the consequence for a
  # likelihood-ratio test is identical.
  collapsed <- function(fit)
    !is.null(fit$degenerate) || length(fit$growth$boundary) > 0L
  degenerate <- c(if (collapsed(restricted)) "restricted",
                  if (collapsed(full))       "full")
  if (length(degenerate))
    warning(sprintf(
      paste0("The %s model carries a collapsed class variance, so this test ",
             "is not interpretable in either direction: a degenerate fit's ",
             "log-likelihood reflects a spike in the likelihood rather than ",
             "fit to the data. Refit it cleanly and test again -- the fit's ",
             "own warning names the remedies. Whichever you choose, apply it ",
             "to both models: a likelihood-ratio test between two fits ",
             "estimated under different constraints is not a test of ",
             "anything."),
      paste(degenerate, collapse = " and ")), call. = FALSE)

  if (b$ll < a$ll)
    warning("`full` has a lower log-likelihood than `restricted`, even ",
            "though it nests it, so the statistic below is negative and is ",
            "not a valid test. It means one of two things. Either `full`'s ",
            "optimizer missed the better solution `restricted` already found, ",
            "in which case refit `full` with a larger n_init; or `restricted` ",
            "has converged on a degenerate solution that outscores the full ",
            "model without describing the data better, in which case inspect ",
            "the fitted variances of both models for a class variance close ",
            "to zero.", call. = FALSE)

  stat <- -2 * (a$ll - b$ll)
  df   <- b$n_params - a$n_params
  out <- list(
    statistic = stat, df = df,
    p_value = if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_,
    ll_restricted = a$ll, ll_full = b$ll,
    params_restricted = a$n_params, params_full = b$n_params
  )
  class(out) <- "lr_test"
  out
}

#' Likelihood-Ratio Test for Two Nested Models (deprecated name)
#'
#' @description
#' Deprecated. The test was never specific to longitudinal models — it accepts
#' any pair of nested fits, and most uses of it in this package are
#' cross-sectional multiple-group invariance tests. Use [`lr_test()`].
#'
#' @param restricted The more constrained model (fewer parameters).
#' @param full The less constrained model.
#' @return A list of class `"lr_test"`.
#' @seealso [`lr_test()`]
#' @export
longitudinal_lrt <- function(restricted, full) {
  .Deprecated("lr_test")
  lr_test(restricted, full)
}

#' @export
print.lr_test <- function(x, ...) {
  cat("\nLikelihood-ratio test for nested models\n")
  cat("---------------------------------------------------------\n")
  cat(sprintf("  Restricted : LL = %12.4f   parameters = %d\n",
              x$ll_restricted, x$params_restricted))
  cat(sprintf("  Full       : LL = %12.4f   parameters = %d\n",
              x$ll_full, x$params_full))
  cat(sprintf("  -2 x diff  : %.4f   df = %d   p = %s\n",
              x$statistic, x$df,
              format.pval(x$p_value, digits = 4, eps = 1e-16)))
  if (!is.na(x$p_value) && x$p_value < 0.05)
    cat("  The restriction is rejected: the full model fits significantly better.\n")
  else if (!is.na(x$p_value))
    cat("  The restriction is not rejected; prefer the more parsimonious model.\n")
  cat("\n")
  invisible(x)
}

# ------------------------------------------------------------------------------
# Display
# ------------------------------------------------------------------------------

#' @export
print.lta_model <- function(x, ...) {
  lg <- x$longitudinal
  K  <- x$n_statuses
  cat("\n")
  cat("=========================================================\n")
  cat("             LATENT TRANSITION ANALYSIS\n")
  cat("=========================================================\n")
  cat(sprintf("Latent statuses    : %d\n", K))
  if ((x$n_classes %||% 1L) > 1L)
    cat(sprintf("Latent classes     : %d%s\n", x$n_classes,
                if (isTRUE(x$mover_stayer))
                  ", the last restricted to no change (mover-stayer)" else ""))
  cat(sprintf("Items x Occasions  : %d x %d\n", lg$n_items, lg$n_times))
  cat(sprintf("Item parameters    : %s, %s\n",
              if (is.list(lg$measurement)) "mixed" else lg$measurement,
              switch(lg$time_invariance,
                     full    = "held equal across occasions",
                     none    = "estimated separately at each occasion",
                     partial = sprintf("%d of %d held equal across occasions",
                                       length(lg$invariant_items), lg$n_items))))
  cat(sprintf("Transitions        : %s\n",
              if (isTRUE(x$tau_homogeneous))
                "one table, shared by every pair of occasions"
              else sprintf("%d tables, one per pair of occasions",
                           lg$n_times - 1L)))
  if (!is.null(x$group_info))
    cat(sprintf("Grouping variable  : %d groups (%s), affecting %s\n",
                length(x$group_info$levels),
                paste(x$group_info$levels, collapse = ", "),
                switch(x$group_effects,
                       both = "prevalences and transitions",
                       initial = "prevalences only",
                       transitions = "transitions only",
                       none = "nothing (fully constrained)")))
  if (!is.null(x$Z_delta) || !is.null(x$Z_tau))
    cat(sprintf("Covariates         : %s initial status, %s transitions\n",
                if (is.null(x$Z_delta)) "no" else
                  as.character(ncol(x$Z_delta) - 1L),
                if (is.null(x$Z_tau)) "no" else
                  as.character(ncol(x$Z_tau) - 1L)))
  cat(sprintf("Converged          : %s (in %d iterations)\n",
              x$converged, x$n_iter))
  # Deleted cases print ahead of the attrition line, which describes only the
  # analysed cases, so the reported n can be reconciled with the input data.
  if (isTRUE(x$missing_data$n_empty_rows > 0L))
    cat(sprintf("Cases Removed      : %d of %d with no observed indicator (n = %d analysed)\n",
                x$missing_data$n_empty_rows, x$missing_data$n_input_rows,
                x$missing_data$n_input_rows - x$missing_data$n_empty_rows))
  if (any(lg$wave_missing))
    cat(sprintf("Wave attrition     : %d case-occasions carried by FIML\n",
                sum(lg$wave_missing)))
  if (isTRUE(x$has_survey_design))
    cat("Survey design      : design-based (linearization) standard errors\n")
  if (identical(x$weight_type, "frequency"))
    cat(sprintf("Case weights       : frequency counts (%s cases in %d rows)\n",
                format(x$n_eff), length(x$weights_vec)))
  cat("---------------------------------------------------------\n")
  .print_fit_indices(
    x$metrics,
    entropy_note = if (is.null(x$metrics$class_entropy)) "" else " (status)")
  if (!is.null(x$metrics$class_entropy))
    cat(sprintf("                   %.4f (class)\n", x$metrics$class_entropy))
  # The quiet channel the mixture models have had all along: how many restarts
  # found this solution. It sits with the other indented metrics rather than in
  # the header block above, which is left-aligned.
  .print_replication_note(x)
  # The boundary from the other side. The note below already names the cells the
  # data have driven to zero; this names the row the prior is holding off it,
  # and only when the prior is carrying enough of that row to matter.
  sm <- .lta_worst_smoothing_row(x$smoothing_influence)
  if (!is.null(sm) && sm$pull > .lta_smoothing_tol)
    cat(sprintf("  Max prior pull : %.2f on %s (%.1f cases expected)\n",
                sm$pull, sm$where, sm$n_expected))
  cat("---------------------------------------------------------\n")
  if ((x$n_classes %||% 1L) > 1L) {
    cat("Latent class sizes:\n")
    cw <- x$class_weights
    names(cw) <- .lta_class_labels(x)
    print(round(cw, 4))
    cat("\nLatent status prevalences by occasion (whole sample):\n")
  } else {
    cat("Latent status prevalences by occasion:\n")
  }
  print(round(.lta_bare_prevalences(x$prevalences), 4))
  if (!is.null(x$boundary))
    cat(sprintf("\nNote: %d transition(s) estimated at the zero boundary; ",
                nrow(x$boundary)),
        "see $boundary.\n", sep = "")
  if ((x$n_classes %||% 1L) > 1L && is.null(x$se))
    cat("\nNote: standard errors are not available for a mixture over chains.\n",
        "      Every parameter here is a point estimate.\n", sep = "")
  cat("=========================================================\n")
  cat("Type summary(model) for transitions, measurement_summary(model) for items.\n")
  invisible(x)
}

#' Summarise a Fitted Latent Transition Model
#'
#' @param object An object returned by [`fit_lta()`].
#' @param digits Number of digits to print.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @export
summary.lta_model <- function(object, digits = 3, ...) {
  K  <- object$n_statuses
  C  <- object$n_classes %||% 1L
  st <- paste0("Status ", seq_len(K))

  cat("\n=========================================================\n")
  cat("     LATENT TRANSITION MODEL - STRUCTURAL PARAMETERS\n")
  cat("=========================================================\n\n")

  if (C > 1L) {
    cat("LATENT CLASS SIZES (pi)\n")
    print(data.frame(Class = .lta_class_labels(object),
                     Size = round(object$class_weights, digits)),
          row.names = FALSE)
    cat("\n")
  }

  cat("INITIAL LATENT STATUS PREVALENCES (delta)\n")
  if (C > 1L) {
    d <- as.data.frame(round(object$delta, digits))
    rownames(d) <- .lta_class_labels(object)
    print(d)
  } else {
    d <- data.frame(Status = st, Prevalence = round(object$delta, digits))
    if (!is.null(object$se$prob_se$delta))
      d$SE <- round(object$se$prob_se$delta, digits)
    print(d, row.names = FALSE)
  }

  cat("\nTRANSITION PROBABILITIES (tau)\n")
  cat("Rows: status at the earlier occasion. Columns: status at the later one.\n")
  if (isTRUE(object$tau_homogeneous))
    cat("A single matrix is shared by every pair of adjacent occasions.\n")

  show_tables <- function(tm, indent = "  ") {
    if (is.matrix(tm)) tm <- list(`all occasions` = tm)
    for (nm in names(tm)) {
      cat(sprintf("\n%s%s\n", indent, nm))
      print(round(tm[[nm]], digits))
    }
  }
  if (C > 1L) {
    labs <- .lta_class_labels(object)
    for (c in seq_len(C)) {
      cat(sprintf("\n%s\n", labs[c]))
      show_tables(transition_matrix(object, class = c), indent = "    ")
    }
  } else {
    show_tables(transition_matrix(object))
  }

  if (!is.null(object$se$prob_se)) {
    tau_se <- object$se$prob_se[grepl("^tau", names(object$se$prob_se))]
    if (length(tau_se)) {
      cat("\n  Standard errors of the transition probabilities\n")
      se_mat <- do.call(rbind, tau_se)
      dimnames(se_mat) <- list(names(tau_se), st)
      print(round(se_mat, digits))
      if (isTRUE(object$se$design_based))
        cat("  Design-based (linearization) standard errors.\n")
      if (isTRUE(object$se$conditional))
        cat("  Conditional on the measurement parameters.\n")
    }
  }

  if (!is.null(object$boundary)) {
    cat("\nBOUNDARY TRANSITIONS (estimated at zero)\n")
    b <- object$boundary
    for (i in seq_len(nrow(b)))
      cat(sprintf("  %soccasion %d: status %d -> status %d\n",
                  if (is.null(b$class)) "" else sprintf("class %d, ", b$class[i]),
                  b$occasion[i], b$from[i], b$to[i]))
    cat("  These cells still consume a degree of freedom in the parameter count.\n")
  }

  cat("\n=========================================================\n")
  invisible(object)
}

#' @export
measurement_summary.lta_model <- function(object, ...) {
  lg  <- object$longitudinal
  inv <- lg$time_invariance == "full"
  cat("\n=========================================================\n")
  cat("        MEASUREMENT MODEL PARAMETERS (by occasion)\n")
  cat("=========================================================\n")
  if (inv)
    cat("Held equal across occasions; one set shown.\n")
  occasions <- if (inv) 1L else seq_len(lg$n_times)
  for (t in occasions) {
    if (!inv) cat(sprintf("\n--- %s ---\n", lg$time_labels[t]))
    sub <- object$mm$models[[t]]
    tmp <- list(mm = sub, n_components = object$n_statuses,
                missing_data = list(any_missing = FALSE))
    # The indicators this occasion's parameters were fitted on, so the table can
    # show the observed marginal beside them. The slice has to be taken and
    # renamed rather than the whole wide matrix passed: the sub-models all carry
    # the *first* occasion's item names, so matching by name against the wide
    # data would report occasion 1's marginals under every occasion's heading.
    slice <- .lta_occasion_indicators(object, lg, sub, t, invariant = inv)
    if (!is.null(slice)) {
      tmp$data           <- slice$data
      tmp$sample_weights <- slice$weights
    }
    class(tmp) <- "mixture_model"
    measurement_summary(tmp)
  }
  invisible(object)
}

# One occasion's columns of the wide indicator matrix, renamed to the names its
# sub-model's parameters carry. Under full time invariance the parameters were
# estimated from every occasion at once, so the marginal they should be read
# against is the one over all of them, and the occasions are stacked rather than
# sliced - one long column per item, with the case weights repeated to match.
#
# The wide layout is time-major and is normalised on input (see
# R/longitudinal_data.R), so occasion t owns columns (t-1)*J + 1 ... t*J. The
# guard on the total width is what makes that arithmetic safe to rely on here;
# anything unexpected returns NULL and the marginal is simply not shown.
.lta_occasion_indicators <- function(object, lg, sub, t, invariant = FALSE) {
  X <- object$data
  if (is.null(X)) return(NULL)
  X <- as.matrix(X)

  pars <- sub$parameters$pis %||% sub$parameters$means
  if (is.null(pars)) return(NULL)
  nm <- colnames(pars)
  J  <- if (!is.null(sub$max_val)) ncol(pars) / sub$max_val else ncol(pars)
  if (J != round(J) || is.null(nm)) return(NULL)

  n_times <- lg$n_times
  if (ncol(X) != n_times * J) return(NULL)

  w  <- object$sample_weights %||% rep(1, nrow(X))
  ts <- if (invariant) seq_len(n_times) else t

  parts <- lapply(ts, function(tt)
    X[, ((tt - 1L) * J + 1L):(tt * J), drop = FALSE])
  out <- do.call(rbind, parts)
  colnames(out) <- if (!is.null(sub$max_val)) sub$item_names else nm
  if (is.null(colnames(out))) return(NULL)

  list(data = out, weights = rep(w, times = length(ts)))
}

#' Plots for a Fitted Latent Transition Model
#'
#' @description
#' Three views, all in base graphics:
#' \describe{
#'   \item{`"prevalence"`}{latent status prevalence across occasions - the
#'     summary of where people are over time;}
#'   \item{`"transitions"`}{a shaded matrix of transition probabilities, one
#'     panel per pair of adjacent occasions, with the values printed in;}
#'   \item{`"profiles"`}{item-response probabilities (or means) by status, which
#'     is what the status labels rest on.}
#' }
#'
#' With more than one latent class, `"prevalence"` draws one panel per class -
#' the classic longitudinal profile plot, in which a stayer class is the
#' flat one - and `"transitions"` one panel per class per pair of occasions.
#' `class` restricts either to a single class.
#'
#' @param x An object returned by [`fit_lta()`].
#' @param type Which view to draw.
#' @param main Plot title.
#' @param status_labels Optional labels for the latent statuses.
#' @param colors Optional colour vector, recycled across statuses.
#' @param class Optional latent class, for a model fitted with `n_classes` > 1
#'   or `mover_stayer = TRUE`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @importFrom graphics par matplot axis legend mtext image box text plot.new
#' @export
plot.lta_model <- function(x, type = c("prevalence", "transitions", "profiles"),
                           main = NULL, status_labels = NULL, colors = NULL,
                           class = NULL, ...) {
  type <- match.arg(type)
  K    <- x$n_statuses
  C    <- x$n_classes %||% 1L
  lg   <- x$longitudinal
  cols <- if (is.null(colors)) rep(.okabe_ito, length.out = K)
          else rep(colors, length.out = K)
  labs <- if (is.null(status_labels)) paste("Status", seq_len(K))
          else status_labels
  if (!is.null(class) && C == 1L)
    stop("`class` applies only to a model with more than one latent class.",
         call. = FALSE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (type == "prevalence") {
    # One panel per class, since a mixture's whole point is that the classes
    # move differently; with one class this is the single panel it always was.
    panels <- if (C == 1L) list(x$prevalences)
              else if (!is.null(class))
                attr(x$prevalences, "by_class")[class]
              else attr(x$prevalences, "by_class")
    ptitles <- if (C == 1L) "" else names(panels) %||% ""
    if (C > 1L && is.null(class)) ptitles <- .lta_class_labels(x)

    # Legend in a reserved strip under the axis: placing it beside the plot at a
    # fixed offset in data units puts it half a panel away when there are two
    # occasions and off the device when there are twenty.
    leg_cols <- max(1L, min(K, 4L))
    leg_rows <- ceiling(K / leg_cols)
    np <- length(panels)
    nr <- if (np == 1L) 1L else ceiling(sqrt(np))
    nc <- ceiling(np / nr)
    par(mfrow = c(nr, nc), mar = c(5, 4, 4, 2),
        oma = c(1.4 * leg_rows, 0, if (np > 1L) 2.5 else 0, 0))
    for (i in seq_along(panels)) {
      matplot(seq_len(lg$n_times), panels[[i]], type = "b",
              pch = rep(15:20, length.out = K),
              lty = 1, lwd = 2, col = cols, ylim = c(0, 1), xaxt = "n",
              xlab = "Occasion", ylab = "Prevalence",
              main = if (np > 1L) ptitles[i] else
                main %||% "Latent status prevalence over time",
              bty = "l", las = 1)
      axis(1, at = seq_len(lg$n_times), labels = lg$time_labels)
    }
    if (np > 1L)
      mtext(main %||% "Latent status prevalence over time, by class",
            outer = TRUE, cex = 1.1, font = 2)
    # Drawn over the whole device rather than at an offset from the plot region,
    # so it lands in the reserved strip whatever the margins work out to be.
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
        mfrow = c(1, 1), new = TRUE)
    plot.new()
    legend("bottom", legend = labs, col = cols,
           pch = rep(15:20, length.out = K), lty = 1, lwd = 2, bty = "n",
           ncol = leg_cols, cex = 0.9, xpd = TRUE)
    return(invisible(x))
  }

  if (type == "transitions") {
    tm <- if (C == 1L) transition_matrix(x)
          else if (!is.null(class)) transition_matrix(x, class = class)
          else {
            per <- transition_matrix(x)
            stats::setNames(
              unlist(lapply(per, function(z)
                if (is.matrix(z)) list(z) else z), recursive = FALSE),
              unlist(lapply(seq_along(per), function(i) {
                z <- per[[i]]
                nm <- if (is.matrix(z)) "all occasions" else names(z)
                paste0(names(per)[i], ": ", nm)
              })))
          }
    if (is.matrix(tm)) tm <- list(`Shared transition matrix` = tm)
    n_panel <- length(tm)
    nr <- ceiling(sqrt(n_panel)); nc <- ceiling(n_panel / nr)
    par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
    shades <- grDevices::colorRampPalette(c("#FFFFFF", "#0072B2"))(64)
    for (nm in names(tm)) {
      m <- tm[[nm]]
      # image() draws the first row at the bottom; reverse so that origin
      # status 1 appears in the top-left, as printed.
      image(seq_len(K), seq_len(K), t(m[K:1, , drop = FALSE]),
            col = shades, zlim = c(0, 1), axes = FALSE,
            xlab = "To", ylab = "From", main = nm)
      axis(1, at = seq_len(K), labels = seq_len(K), tick = FALSE)
      axis(2, at = seq_len(K), labels = rev(seq_len(K)), tick = FALSE, las = 1)
      box()
      for (i in seq_len(K)) for (j in seq_len(K))
        text(j, K - i + 1, sprintf("%.2f", m[i, j]),
             col = if (m[i, j] > 0.5) "white" else "grey20", cex = 0.9)
    }
    mtext(main %||% "Transition probabilities", outer = TRUE, line = 1,
          cex = 1.1, font = 2)
    return(invisible(x))
  }

  # type == "profiles"
  # Time 1 status prevalences rather than `delta`, which gains a class index in
  # a mixture; for a single chain the two are the same vector.
  tmp <- list(mm = x$mm$models[[1]], n_components = K,
              weights = .lta_bare_prevalences(x$prevalences)[1, ], data = NULL)
  class(tmp) <- "mixture_model"
  plot.mixture_model(tmp,
                     main = main %||% "Item responses by latent status",
                     class_labels = labs, colors = cols)
  invisible(x)
}

#' @export
classification_diagnostics.lta_model <- function(object, ...) {
  K <- object$n_statuses
  w <- object$weights_vec
  cat("\n=========================================================\n")
  cat("     AVERAGE POSTERIOR PROBABILITIES, BY OCCASION\n")
  cat("=========================================================\n")
  cat("Rows: modal status assignment | Columns: mean probability\n")
  if ((object$n_classes %||% 1L) > 1L)
    cat("Statuses only: how cleanly the latent *classes* separate is not\n",
        "reported here; see the class entropy in print().\n", sep = "")
  out <- vector("list", object$longitudinal$n_times)
  for (t in seq_len(object$longitudinal$n_times)) {
    g <- object$gamma[[t]]
    m <- .ave_pp(g, w, K)
    dimnames(m) <- list(paste("Assigned", seq_len(K)),
                        paste("Prob", seq_len(K)))
    cat(sprintf("\n--- %s ---\n", object$longitudinal$time_labels[t]))
    print(round(m, 3))
    out[[t]] <- list(ave_pp = m, table = .classification_table(g, w, K))
  }
  cat("=========================================================\n")

  # A status assignment is made at every occasion, so the classification error
  # is reported per occasion rather than pooled: the measurement model may be
  # invariant, but the status distribution it is applied to is not.
  cat("\nClassification error by occasion:\n")
  for (t in seq_along(out))
    cat(sprintf("  %-12s %.4f\n", object$longitudinal$time_labels[t],
                attr(out[[t]]$table, "error")))
  cat("=========================================================\n")

  names(out) <- object$longitudinal$time_labels
  invisible(out)
}
