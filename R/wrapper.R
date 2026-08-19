# ==============================================================================
# S3 User Wrappers and Pipeline Tools (LCA/LPA Mixture Engine)
# ==============================================================================

sort_model_classes <- function(model_state) {
  K <- model_state$n_components
  if (K <= 1) return(model_state)

  new_order <- order(model_state$weights, decreasing = TRUE)
  model_state$weights <- model_state$weights[new_order]

  # --- Sort flat measurement model parameters ---
  if (!is.null(model_state$mm$parameters[["pis"]])) {
    model_state$mm$parameters$pis <-
      model_state$mm$parameters$pis[new_order, , drop = FALSE]
  }
  # Gaussian LPA: also sort means and covariances
  if (!is.null(model_state$mm$parameters[["means"]])) {
    model_state$mm$parameters$means <-
      model_state$mm$parameters$means[new_order, , drop = FALSE]
  }
  if (!is.null(model_state$mm$parameters[["covariances"]])) {
    model_state$mm$parameters$covariances <-
      model_state$mm$parameters$covariances[new_order, , drop = FALSE]
  }
  # Poisson LCA: class-specific rates are indexed by class like pis and means
  if (!is.null(model_state$mm$parameters[["rates"]])) {
    model_state$mm$parameters$rates <-
      model_state$mm$parameters$rates[new_order, , drop = FALSE]
  }
  # LCGA: growth coefficients are one row per class, and the gaussian family's
  # residual variance is one entry per class.
  if (!is.null(model_state$mm$parameters[["coefs"]])) {
    model_state$mm$parameters$coefs <-
      model_state$mm$parameters$coefs[new_order, , drop = FALSE]
    model_state$mm$parameters$dispersion <-
      model_state$mm$parameters$dispersion[new_order]
  }
  # GMM: growth-factor means are one row per class, residual variances one row
  # per class, and the growth-factor covariance one matrix per class. All are
  # stored per class even where a constraint makes the classes share a value, so
  # all are permuted unconditionally. The growth-factor regressions on
  # covariates are stored the same way and join them when there are any.
  if (!is.null(model_state$mm$parameters[["alpha"]])) {
    model_state$mm$parameters$alpha <-
      model_state$mm$parameters$alpha[new_order, , drop = FALSE]
    model_state$mm$parameters$theta <-
      model_state$mm$parameters$theta[new_order, , drop = FALSE]
    model_state$mm$parameters$psi <-
      model_state$mm$parameters$psi[new_order]
    if (!is.null(model_state$mm$parameters[["gamma"]]))
      model_state$mm$parameters$gamma <-
        model_state$mm$parameters$gamma[new_order]
  }

  # --- Sort nested measurement model sub-model parameters ---
  # The flat-parameter block above only touches model_state$mm$parameters, which
  # is empty for nested models.  Sub-model parameters live one level deeper at
  # model_state$mm$models[[name]]$parameters and must be sorted independently.
  if (inherits(model_state$mm, c("nested", "blocks"))) {
    for (name in names(model_state$mm$models)) {
      sub <- model_state$mm$models[[name]]
      if (!is.null(sub$parameters[["pis"]]))
        sub$parameters$pis <- sub$parameters$pis[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["means"]]))
        sub$parameters$means <- sub$parameters$means[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["covariances"]]))
        sub$parameters$covariances <-
          sub$parameters$covariances[new_order, , drop = FALSE]
      if (!is.null(sub$parameters[["rates"]]))
        sub$parameters$rates <- sub$parameters$rates[new_order, , drop = FALSE]
      model_state$mm$models[[name]] <- sub
    }
  }

  # --- Sort structural model parameters ---
  if (!is.null(model_state$sm)) {

    sort_sm_params <- function(sm) {
      if (!is.null(sm$parameters[["beta"]])) {
        sm$parameters$beta <- sm$parameters$beta[new_order, , drop = FALSE]
        if (!is.null(sm$parameters[["hessian"]])) {
          H       <- sm$parameters$hessian
          D       <- ncol(sm$parameters$beta)
          idx_map <- as.vector(sapply(new_order, function(k) ((k-1)*D + 1):(k*D)))
          sm$parameters$hessian <- H[idx_map, idx_map, drop = FALSE]
        }
        # The survey-robust covariance is blocked by class in the same layout
        # as the Hessian, so it is permuted with the identical index map.
        if (!is.null(sm$parameters[["V_robust"]])) {
          Vr      <- sm$parameters$V_robust
          D       <- ncol(sm$parameters$beta)
          idx_map <- as.vector(sapply(new_order, function(k) ((k-1)*D + 1):(k*D)))
          sm$parameters$V_robust <- Vr[idx_map, idx_map, drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["pis"]])) {
        sm$parameters$pis <- sm$parameters$pis[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["means"]])) {
        sm$parameters$means <- sm$parameters$means[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["covariances"]])) {
        sm$parameters$covariances <-
          sm$parameters$covariances[new_order, , drop = FALSE]
      }
      if (!is.null(sm$parameters[["ses"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$ses[1, 1:K] <- sm$parameters$ses[1, new_order]
        } else {
          sm$parameters$ses <- sm$parameters$ses[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["Sigma_mu"]])) {
        sm$parameters$Sigma_mu <-
          sm$parameters$Sigma_mu[new_order, new_order, drop = FALSE]
      }
      if (!is.null(sm$parameters[["cov_theta"]])) {
        # cov_theta is L x L where the first K rows/cols are intercepts.
        # Reorder only the intercept block; slope block stays unchanged.
        K_ct  <- sm$n_components
        L_ct  <- nrow(sm$parameters$cov_theta)
        D_ct  <- L_ct - K_ct
        idx   <- c(new_order, if (D_ct > 0) (K_ct + seq_len(D_ct)) else integer(0))
        sm$parameters$cov_theta <-
          sm$parameters$cov_theta[idx, idx, drop = FALSE]
      }
      # Reorder the robust-sandwich hessian for distal_pooled / distal_regression.
      # This hessian is K x K (or larger) and is indexed by class in the
      # pre-sort order.  After sort_model_classes reorders beta_pooled, the
      # hessian must be permuted to match so that Sigma = pinv(-H) is aligned.
      if (!is.null(sm$parameters[["hessian"]]) &&
          inherits(sm, c("distal_pooled", "distal_regression"))) {
        H_sm  <- sm$parameters$hessian
        K_sm  <- sm$n_components
        L_sm  <- nrow(H_sm)
        D_sm  <- L_sm - K_sm
        idx_h <- c(new_order,
                   if (D_sm > 0) (K_sm + seq_len(D_sm)) else integer(0))
        sm$parameters$hessian <- H_sm[idx_h, idx_h, drop = FALSE]
      }
      if (!is.null(sm$parameters[["betas"]])) {
        if (length(dim(sm$parameters$betas)) == 3) {
          sm$parameters$betas <- sm$parameters$betas[new_order, , , drop = FALSE]
        } else {
          sm$parameters$betas <- sm$parameters$betas[new_order, , drop = FALSE]
        }
      }
      if (!is.null(sm$parameters[["beta_pooled"]])) {
        if (inherits(sm, "distal_continuous_pooled")) {
          sm$parameters$beta_pooled[1, 1:K] <- sm$parameters$beta_pooled[1, new_order]
        } else {
          # Guard: when beta_pooled is a degenerate 0x0 placeholder (produced by
          # a constant-outcome distal_pooled model), there is nothing to reorder.
          if (nrow(sm$parameters$beta_pooled) > 0 &&
              ncol(sm$parameters$beta_pooled) >= K) {
            sm$parameters$beta_pooled[, 1:K] <-
              sm$parameters$beta_pooled[, new_order, drop = FALSE]
          }
        }
      }
      return(sm)
    }

    if (inherits(model_state$sm, "nested")) {
      for (name in names(model_state$sm$models))
        model_state$sm$models[[name]] <- sort_sm_params(model_state$sm$models[[name]])
    } else {
      model_state$sm <- sort_sm_params(model_state$sm)
    }
  }

  if (!is.null(model_state$log_resp))
    model_state$log_resp <- model_state$log_resp[, new_order, drop = FALSE]

  return(model_state)
}

#' Print Measurement Model Parameters
#'
#' @description
#' Prints a formatted table of the fitted measurement model parameters:
#' item-response probabilities for categorical models, or means for Gaussian
#' models. Results are broken down by latent class. Handles both flat and
#' nested (mixed) measurement models.
#'
#' For a growth model — \code{\link{fit_gmm}} or \code{\link{fit_lcga}} — the
#' measurement parameters are the growth-factor means, their variances and
#' covariances, the residual variances and the fitted trajectory, and those are
#' what the table holds. A parameter held equal across classes is repeated once
#' per class rather than reported once, so the table can be joined to anything
#' else indexed by class; the constraint is stated in the printed heading.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param scale For categorical indicators, what scale the item parameters are
#'   reported on. \code{"probability"} (the default) is unchanged from before
#'   this argument existed. \code{"logit"} reports \code{qlogis()} of the same
#'   table. \code{"effect"} reports the effect-coded parameterisation several
#'   other programs use by default -- an item intercept plus one deviation per
#'   class, the deviations summing to zero -- which is what lets a mixtureEM
#'   measurement model be placed beside such a program's printed output;
#'   binary indicators only, since a polytomous item's effect coding is a
#'   modelling choice (ordinal with fixed scores, giving one class effect per
#'   class, versus nominal, giving one per category) that the package does
#'   not make for you, and a polytomous item under \code{"effect"} is refused
#'   with an error rather than guessed at. The \code{overall} column is
#'   dropped on the \code{"logit"} and \code{"effect"} scales, since the
#'   observed marginal has no meaningful transform there. Ignored for
#'   continuous means and count rates, which have no probability to rescale.
#' @param ... Passed to methods.
#'
#' @return Invisibly, a data frame in long format with one row per item,
#'   response category (polytomous items only, \code{NA} otherwise), and
#'   class: columns \code{block} (sub-model name for mixed measurement
#'   models, \code{NA} otherwise), \code{parameter} (\code{"probability"},
#'   \code{"mean"}, or \code{"rate"}; for a growth model
#'   \code{"growth_mean"}, \code{"growth_variance"},
#'   \code{"growth_covariance"}, \code{"growth_regression"},
#'   \code{"residual_variance"} or \code{"fitted"}), \code{item},
#'   \code{category}, \code{class}, \code{estimate}, and \code{overall}. The
#'   same numbers are printed as formatted tables.
#'
#'   \code{overall} is the observed marginal for the item — the weighted sample
#'   proportion beside a probability, the weighted sample mean beside a mean or
#'   a rate — repeated down the class rows so the frame stays joinable on
#'   \code{class}. It is what makes a conditional number readable: a class
#'   endorsing an item at .62 is unremarkable when the sample sits at .60 and is
#'   most of what defines the class when the sample sits at .12. It is
#'   \code{NA}, and the printed column is dropped, for a fit that does not store
#'   its raw indicators or whose item parameters cannot be matched to them by
#'   name.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' measurement_summary(fit)
#' params <- measurement_summary(fit)   # reuse the table programmatically
#'
#' @export
measurement_summary <- function(object, ...) UseMethod("measurement_summary")

#' @rdname measurement_summary
#' @export
measurement_summary.default <- function(object,
                                        scale = c("probability", "logit",
                                                  "effect"),
                                        ...) {
  scale <- match.arg(scale)
  if (identical(scale, "effect")) {
    poly <- .polytomous_block_names(object$mm)
    if (length(poly))
      stop("`scale = \"effect\"` is not supported yet for a polytomous ",
           "item: ", paste(poly, collapse = ", "), ". A polytomous item's ",
           "effect coding depends on whether it is ordinal (one class ",
           "effect per class) or nominal (one per category), which is a ",
           "modelling choice this function does not make for you. Use ",
           "`scale = \"logit\"` instead, or report probabilities.",
           call. = FALSE)
  }
  K <- object$n_components
  cat("=========================================================\n")
  cat("             MEASUREMENT MODEL PARAMETERS                \n")
  cat("=========================================================\n")

  collected      <- list()
  # Set by any block whose sample marginal could not be worked out, so the note
  # under the table is printed once however many blocks the model has.
  any_no_overall <- FALSE

  print_item_matrix <- function(mat, title, sub_model = NULL,
                                parameter = "estimate",
                                block = NA_character_,
                                intercept = NULL, show_overall = TRUE) {
    cat(sprintf("\n%s\n", title))
    item_names <- colnames(mat)
    base_items <- NULL
    categories <- NULL

    if (is.null(item_names)) {
      if (!is.null(sub_model) && !is.null(sub_model$max_val)) {
        M       <- sub_model$max_val
        n_items <- ncol(mat) / M
        base    <- sub_model$item_names
        if (is.null(base) || length(base) != n_items)
          base <- paste0("Poly_Item_", seq_len(n_items))
        item_names <- paste0(rep(base, each = M),
                             " (Cat ", rep(seq_len(M), times = n_items), ")")
        base_items <- rep(base, each = M)
        categories <- rep(seq_len(M), times = n_items)
      } else {
        item_names <- paste0("Item_", 1:ncol(mat))
      }
    }
    if (is.null(base_items)) base_items <- item_names
    if (is.null(categories)) categories <- rep(NA_integer_, length(item_names))

    # Long indicator names widen the table only up to the shortening cap;
    # anything longer is abbreviated with a key printed under the table. The
    # returned data frame keeps the full names.
    disp    <- .shorten_labels(item_names, width = 30L)
    label_w <- max(20L, max(nchar(disp)))

    # The sample marginal goes between the label and the classes, so each row
    # reads as "here is the item overall, and here is how each class departs
    # from it" - which is the comparison the class parameters are for. Dropped
    # entirely, rather than filled with NAs, when it cannot be computed -- or
    # when the caller has already said it is meaningless on this scale.
    overall <- if (show_overall) .observed_marginals(object, mat, sub_model,
                                                      parameter) else NULL
    if (!is.null(overall) && length(overall) != ncol(mat)) overall <- NULL
    if (show_overall && is.null(overall)) any_no_overall <<- TRUE

    cat(sprintf("%-*s", label_w, "Indicator"))
    if (!is.null(intercept)) cat(" | Intercept")
    if (!is.null(overall)) cat(" | Overall")
    for (k in 1:K) cat(sprintf(" | Class %d", k))
    cat("\n")
    cat(paste0(rep("-", label_w + (K + !is.null(overall) +
                                    !is.null(intercept)) * 10),
               collapse = ""), "\n")

    for (j in 1:ncol(mat)) {
      cat(sprintf("%-*s", label_w, disp[j]))
      if (!is.null(intercept)) cat(sprintf(" | %7.3f", intercept[j]))
      if (!is.null(overall)) cat(sprintf(" | %7.3f", overall[j]))
      for (k in 1:K) cat(sprintf(" | %7.3f", mat[k, j]))
      cat("\n")
    }
    .cat_label_legend(disp, indent = "")

    # Long-format rows for the returned data frame. as.vector(mat) walks the
    # K x J matrix column by column, so class varies fastest within each item.
    J <- ncol(mat)
    collected[[length(collected) + 1L]] <<- data.frame(
      block     = rep(block, K * J),
      parameter = rep(parameter, K * J),
      item      = rep(base_items, each = K),
      category  = rep(categories, each = K),
      class     = rep(seq_len(K), times = J),
      estimate  = as.vector(mat),
      # Constant within item and category, repeated down the K class rows, so
      # the frame stays joinable on `class` and a per-class departure from the
      # marginal is a subtraction rather than a reshape. NA for a block whose
      # marginal could not be worked out; the column itself is always present,
      # since a data frame that gains and loses columns by fit is worse to
      # program against than one that carries an NA.
      overall   = rep(overall %||% rep(NA_real_, J), each = K),
      stringsAsFactors = FALSE)
  }

  # Title suffix and column suppression are the only places `scale` touches
  # anything outside the categorical-probability blocks themselves.
  scale_suffix <- switch(scale, probability = "",
                         logit = " (logit scale)", effect = " (effect-coded)")
  show_overall <- identical(scale, "probability")

  mm <- object$mm
  if (inherits(mm, c("nested", "blocks"))) {
    for (name in names(mm$models)) {
      sub_mm <- mm$models[[name]]
      if (!is.null(sub_mm$parameters$pis)) {
        tr <- .scale_categorical_block(sub_mm$parameters$pis, scale)
        print_item_matrix(tr$mat,
                          paste0("Categorical Probabilities: ", toupper(name),
                                 scale_suffix),
                          sub_mm, "probability", name,
                          intercept = tr$intercept, show_overall = show_overall)
      }
      if (!is.null(sub_mm$parameters$means))
        print_item_matrix(sub_mm$parameters$means,
                          paste("Continuous Means:", toupper(name)),
                          sub_mm, "mean", name)
      if (!is.null(sub_mm$parameters$rates))
        print_item_matrix(sub_mm$parameters$rates,
                          paste("Count Rates:", toupper(name)),
                          sub_mm, "rate", name)
    }
  } else {
    if (!is.null(mm$parameters$pis)) {
      tr <- .scale_categorical_block(mm$parameters$pis, scale)
      print_item_matrix(tr$mat, paste0("CATEGORICAL PROBABILITIES", scale_suffix),
                        mm, "probability",
                        intercept = tr$intercept, show_overall = show_overall)
    }
    if (!is.null(mm$parameters$means))
      print_item_matrix(mm$parameters$means, "CONTINUOUS MEANS", mm, "mean")
    if (!is.null(mm$parameters$rates))
      print_item_matrix(mm$parameters$rates, "COUNT RATES (lambda)", mm,
                        "rate")
  }
  if (any_no_overall)
    cat("\nThe Overall column, holding the observed marginal for each item, is ",
        "omitted above: this fit either does not store its raw indicators - in ",
        "which case refitting with the current version enables it - or holds ",
        "item parameters that cannot be matched to them by name, as a ",
        "multiple-group measurement model does.\n", sep = "")
  if (!is.null(object$missing_data) && isTRUE(object$missing_data$any_missing)) {
    md <- object$missing_data
    cat(sprintf("\nMissing data: %d of %d cells (%.1f%%) across %d item%s, handled via %s.\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  .print_recode_note(object$binary_recode)
  .print_boundary_note(do.call(rbind, collected))
  cat("=========================================================\n")
  invisible(do.call(rbind, collected))
}

# Names of any polytomous categorical blocks in a measurement model, or
# character(0) if every categorical block is binary. A block counts as
# polytomous when it stores a `max_val` -- the Bernoulli branch of
# `.categorical_item_probs()` (R/fit_diagnostics.R) has none, since a binary
# item's stored parameter is a single P(y = 1) rather than one probability
# per category.
.polytomous_block_names <- function(mm) {
  if (inherits(mm, c("nested", "blocks"))) {
    keep <- vapply(mm$models, function(sub)
      !is.null(sub$parameters$pis) && !is.null(sub$max_val), logical(1))
    return(names(mm$models)[keep])
  }
  if (!is.null(mm$parameters$pis) && !is.null(mm$max_val))
    return("the categorical block")
  character(0)
}

# Rescale a K x J block of binary item-response probabilities for
# measurement_summary()'s `scale` argument. "probability" is a no-op.
# "logit" is qlogis() elementwise and needs no restriction on the item type.
# "effect" recovers another program's effect coding for a binary item: an
# item intercept -- half the negative mean logit across classes, matching the
# reported category's intercept row of an equivalent two-category coding,
# whose other category's intercept is this value's negative -- plus one class
# deviation per class, the deviations summing to zero by construction. The
# caller has already refused a polytomous block under "effect" before this is
# reached. Verified against another program's printed output on a reference
# fit: matches to three decimals on both the intercept and the four class
# deviations of a boundary-adjacent item.
.scale_categorical_block <- function(mat, scale) {
  if (identical(scale, "probability")) return(list(mat = mat, intercept = NULL))
  logit <- stats::qlogis(mat)
  if (identical(scale, "logit")) return(list(mat = logit, intercept = NULL))
  mu <- colMeans(logit)
  list(mat = sweep(logit, 2, mu, "-"), intercept = -mu / 2)
}

# The observed marginal for each column of a block of item parameters, on the
# same scale as the parameters themselves: the weighted sample proportion for a
# probability, the weighted sample mean for a mean or a rate.
#
# This is what puts a conditional number in context. A class endorsing an item
# at .62 is unremarkable when everybody endorses it at .60 and is most of what
# defines the class when the sample sits at .12, and the table alone cannot tell
# those apart. The comparison is the reason the column is the *observed*
# marginal rather than the model-implied one, sum_k pi_k theta_jk: a
# model-implied column would be another number the model produced, and would
# agree with the class parameters by construction, which is exactly the property
# that makes it useless as a benchmark.
#
# The indicators are read as stored on the fit, which is after any binary
# recode, so a proportion here is of the same level the probability beside it
# is of. Returns NULL when the raw data is unavailable or the columns cannot be
# lined up; the caller then drops the column rather than printing NAs.
.observed_marginals <- function(object, mat, sub_model = NULL,
                                parameter = "estimate") {
  data <- object$data
  if (is.null(data) || !is.matrix(data) && !is.data.frame(data)) return(NULL)
  data <- as.matrix(data)
  if (nrow(data) == 0L) return(NULL)

  w <- object$sample_weights %||% rep(1, nrow(data))
  if (length(w) != nrow(data)) w <- rep(1, nrow(data))

  # Weighted mean over the observed cases of one column.
  wmean <- function(v) {
    ok <- !is.na(v)
    if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
    sum(w[ok] * v[ok]) / sum(w[ok])
  }

  # Polytomous blocks hold M columns per item, laid out item-major, and the
  # marginal for column (j-1)*M + m is the share of cases answering m to item j.
  # data.matrix() codes the categories 1..M, which is what max_val counts.
  if (identical(parameter, "probability") && !is.null(sub_model$max_val)) {
    M       <- sub_model$max_val
    n_items <- ncol(mat) / M
    if (n_items != round(n_items)) return(NULL)
    base    <- sub_model$item_names
    if (is.null(base) || length(base) != n_items) return(NULL)
    matched <- .match_indicator_columns(base, data)
    if (is.null(matched)) return(NULL)

    out <- numeric(ncol(mat))
    for (j in seq_len(n_items)) {
      v <- data[, matched[j]]
      for (m in seq_len(M))
        out[(j - 1L) * M + m] <- wmean(as.numeric(v == m))
    }
    return(out)
  }

  matched <- .match_indicator_columns(colnames(mat), data)
  if (is.null(matched)) return(NULL)
  vapply(matched, function(cj) wmean(as.numeric(data[, cj])), numeric(1))
}

# Say which level each reported probability belongs to, for items that were not
# supplied as 0/1. Without this the table is ambiguous in the one way that
# matters: "0.87" means nothing until you know 0.87 of what.
.print_recode_note <- function(map) {
  if (is.null(map) || !length(map)) return(invisible(NULL))
  pairs <- vapply(map, function(m) sprintf("%s = %s", m$item, m$one),
                  character(1))
  cat("\nProbabilities are of: ", paste(pairs, collapse = ", "),
      " (recoded to 0/1 on input).\n", sep = "")
  invisible(NULL)
}

# Name the probabilities that have run to 0 or 1.
#
# Two things follow for the reader, and neither is obvious from the table. The
# class is defined partly by an item that nobody in it endorses (or that
# everybody does), which is a substantive fact about the data worth checking.
# And the standard error there is not interpretable: at the boundary the
# information matrix loses rank, so confidence intervals and significance tests
# for that parameter mean nothing (Galindo Garre & Vermunt, 2006). The weak
# Dirichlet prior keeps estimates strictly inside the interval, so a value this
# close to the edge means the likelihood was pushing hard against it.
.print_boundary_note <- function(df, tol = 1e-3) {
  if (is.null(df) || !nrow(df) || !"parameter" %in% names(df)) return(invisible(NULL))
  keep <- df$parameter == "probability" & is.finite(df$estimate) &
    (df$estimate <= tol | df$estimate >= 1 - tol)
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) return(invisible(NULL))

  hit <- df[keep, , drop = FALSE]
  lab <- sprintf("%s in class %s", hit$item, hit$class)
  if (!all(is.na(hit$block))) lab <- sprintf("%s (%s)", lab, hit$block)
  show <- utils::head(lab, 5L)
  cat(sprintf(
    paste0("\nAt the boundary: %s%s. These probabilities have run to 0 or 1, ",
           "so the class is defined partly by an item every case in it gives ",
           "the same answer to, and their standard errors are not ",
           "interpretable. There are two ways on: read it substantively, since ",
           "an item every member of a class answers identically is often the ",
           "finding rather than a fault; or, if that parameter needs a ",
           "standard error, refit with a stronger prior than the default of 1 ",
           "- `bayes_constants = list(categorical = 2)` - which holds the ",
           "estimate off the edge at the cost of shrinking it slightly toward ",
           "the item's marginal.\n"),
    paste(show, collapse = "; "),
    if (length(lab) > 5L) sprintf(" and %d more", length(lab) - 5L) else ""))
  invisible(NULL)
}

#' Print Classification Diagnostics
#'
#' @description
#' Prints the two tables that describe how cleanly the model assigns cases.
#'
#' The **Average Posterior Probability (AvePP)** matrix has one row per set of
#' observations modally assigned to a class and one column per class, holding
#' the mean posterior probability. High values on the diagonal, low values off
#' it, indicate well-separated classes.
#'
#' The **classification table** cross-classifies the probabilistic memberships
#' against the modal assignment, and yields the classification error: the
#' proportion of cases the modal rule is expected to place in the wrong class.
#' See [`classification_table()`] for the details, and [`absolute_fit()`] and
#' [`bivariate_residuals()`] for the fit of the model itself rather than the
#' quality of its assignments.
#'
#' The diagonal of the AvePP matrix reads as a percentage of the cases assigned
#' to that class: a value of 0.91 "suggests that 91% of subjects in the assigned
#' class fit that category, while 9% of the subjects in that class do not
#' accurately fit that category" (Fanti & Henrich, 2010, as reported by Lee et
#' al., 2023, p. 653).
#'
#' This function does not compute entropy, which is the other number solutions
#' are judged by. That is `fit$metrics$entropy`, printed as `Rel. Entropy` and
#' carried as the `Entropy` column of [`compare_mixtures()`] and
#' [`compare_longitudinal()`], where it is documented.
#'
#' Both tables use the case weights when the model was fitted with any.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Passed to methods.
#'
#' @return Invisibly, a list with `ave_pp` (the K x K matrix), `table` (the
#'   classification table) and `error` (the classification error). All are also
#'   printed to the console.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' classification_diagnostics(fit)
#'
#' @references
#' Celeux, G., & Soromenho, G. (1996). An entropy criterion for assessing the
#' number of clusters in a mixture model. \emph{Journal of Classification},
#' \emph{13}(2), 195-212. \doi{10.1007/BF01246098}
#'
#' Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction to
#' growth mixture models (GMM). In \emph{International Encyclopedia of
#' Education} (4th ed., Vol. 14, pp. 646-655). Elsevier.
#' \doi{10.1016/B978-0-12-818630-5.10076-4}
#'
#' Nagin, D. S. (2005). \emph{Group-Based Modeling of Development}. Harvard
#' University Press.
#'
#' @seealso [`class_assignments()`] for the per-case assignment these
#'   diagnostics summarise.
#'
#' @export
classification_diagnostics <- function(object, ...)
  UseMethod("classification_diagnostics")

#' @rdname classification_diagnostics
#' @export
classification_diagnostics.default <- function(object, ...) {
  resp <- exp(object$log_resp)
  K    <- object$n_components
  w    <- object$sample_weights %||% rep(1, nrow(resp))

  ave_pp <- .ave_pp(resp, w, K)
  rownames(ave_pp) <- paste("Assigned Class", 1:K)
  colnames(ave_pp) <- paste("Prob C", 1:K)

  cat("=========================================================\n")
  cat("          AVERAGE POSTERIOR PROBABILITIES (AvePP)        \n")
  cat("=========================================================\n")
  cat("Rows: Modal Assignment | Columns: Mean Probability\n\n")
  print(round(ave_pp, 3))
  cat("=========================================================\n")

  cat("\n")
  tab <- .classification_table(resp, w, K)
  print(tab)

  invisible(list(ave_pp = ave_pp, table = tab, error = attr(tab, "error")))
}

#' Class Sizes of a Fitted Mixture Model
#'
#' @description
#' Returns the estimated size of each latent class in the three forms applied
#' papers report: the model's class proportion, the expected number of cases,
#' and the number of cases modally assigned to the class. Case weights are
#' used when the model was fitted with any.
#'
#' @details
#' The smallest class is one of the things readers judge a solution by, and
#' there are two published conventions for it. Lee et al. (2023, p. 654): "If
#' the smallest class contains less than 5\% of the sample and/or the sample
#' size for the smallest class is less than 25, it is recommended that the model
#' only be retained as the optimal model if the researcher can accurately defend
#' what is gained from this small class given the possibility of low power and a
#' lack of statistical precision." Jung and Wickrama (2008, p. 312) give a
#' weaker floor among their checks: no less than 1\% of the total in any class.
#'
#' Both are reporting conventions, and mixtureEM enforces neither. This function
#' applies no threshold, raises no warning and filters nothing; a small class is
#' estimated and returned like any other. What the conventions ask of you is a
#' defence, not a deletion: keep the class if it can be justified substantively,
#' drop it if it cannot, and report the number either way.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Passed to methods.
#'
#' @seealso [`class_assignments()`] for the per-case assignment the
#'   \code{n_modal} column counts.
#'
#' @return A data frame with one row per class: \code{class},
#'   \code{proportion} (model-estimated class weight), \code{n_expected}
#'   (proportion times the analysed sample size), and \code{n_modal} (cases
#'   assigned by highest posterior probability).
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#' class_sizes(fit)
#'
#' @references
#' Jung, T., & Wickrama, K. A. S. (2008). An introduction to latent class growth
#' analysis and growth mixture modeling. \emph{Social and Personality Psychology
#' Compass}, \emph{2}(1), 302-317. \doi{10.1111/j.1751-9004.2007.00054.x}
#'
#' Lee, T. K., Wickrama, K. A. S., & O'Neal, C. W. (2023). An introduction to
#' growth mixture models (GMM). In \emph{International Encyclopedia of
#' Education} (4th ed., Vol. 14, pp. 646-655). Elsevier.
#' \doi{10.1016/B978-0-12-818630-5.10076-4}
#'
#' @export
class_sizes <- function(object, ...) UseMethod("class_sizes")

#' @rdname class_sizes
#' @export
class_sizes.mixture_model <- function(object, ...) {
  K     <- object$n_components
  resp  <- exp(object$log_resp)
  w     <- object$sample_weights %||% rep(1, nrow(resp))
  modal <- max.col(resp, ties.method = "first")
  n_tot <- sum(w)

  data.frame(
    class      = seq_len(K),
    proportion = as.vector(object$weights),
    n_expected = as.vector(object$weights) * n_tot,
    n_modal    = vapply(seq_len(K), function(k) sum(w[modal == k]), numeric(1))
  )
}

#' Class Assignments for Each Case
#'
#' @description
#' The class each case is assigned to, and the posterior probabilities behind
#' that assignment. This is the accessor for the per-case classification, so
#' that reaching it does not mean indexing into the fitted object's internals.
#'
#' @details
#' Modal class assignment discards classification error. Do not use the returned
#' class as though it were an observed variable in a subsequent regression,
#' ANOVA or t-test: doing so attenuates the association, severely when the
#' classes are not well separated (Bolck, Croon & Hagenaars, 2004; Vermunt,
#' 2010; Bakk, Tekle & Vermunt, 2013). Use [`add_covariates()`] and
#' [`add_outcome()`], which correct for it. This function is for plotting,
#' exporting and describing a solution.
#'
#' @param object A fitted model: a \code{mixture_model} (including the growth
#'   models) or an \code{lta_model}.
#' @param type What to return. \code{"modal"} (default) gives the assigned
#'   class; \code{"posterior"} the full matrix of posterior probabilities;
#'   \code{"both"} a data frame carrying the assignment, its probability, and
#'   the posterior columns.
#' @param ... Passed to methods.
#'
#' @return For \code{"modal"}, an integer vector of length n. For
#'   \code{"posterior"}, an n-by-K matrix with the class labels as column names.
#'   For \code{"both"}, a data frame with \code{class}, \code{probability} (the
#'   assigned class's posterior probability, i.e. a per-case classification
#'   certainty), and then the K posterior columns.
#'
#' @references
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent structure
#' models with categorical variables: One-step versus three-step estimators.
#' \emph{Political Analysis}, \emph{12}(1), 3–27. \doi{10.1093/pan/mph001}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis}, \emph{18}(4), 450–469.
#' \doi{10.1093/pan/mpq025}
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#' between latent class membership and external variables using bias-adjusted
#' three-step approaches. \emph{Sociological Methodology}, \emph{43}(1),
#' 272–311. \doi{10.1177/0081175012470644}
#'
#' @seealso [`class_sizes()`], [`classification_table()`],
#'   [`classification_diagnostics()`].
#'
#' @examples
#' set.seed(1)
#' X   <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#' table(class_assignments(fit))
#' head(class_assignments(fit, "both"))
#' # To relate the classes to an external variable, do not regress on the
#' # assigned class - use the bias-adjusted third step instead:
#' # add_outcome(fit, y)
#'
#' @export
class_assignments <- function(object, type = c("modal", "posterior", "both"),
                              ...)
  UseMethod("class_assignments")

#' @rdname class_assignments
#' @export
class_assignments.default <- function(object,
                                      type = c("modal", "posterior", "both"),
                                      ...) {
  type <- match.arg(type)
  if (is.null(object$log_resp))
    stop("`object` carries no posterior class probabilities, so there is ",
         "nothing to assign from.", call. = FALSE)

  resp <- exp(object$log_resp)
  # ties.method = "first" matches get_modal_resp() in R/corrections.R. The two
  # must not disagree: one is what the user sees, the other is what the modal
  # three-step correction is a table of.
  modal <- max.col(resp, ties.method = "first")
  if (type == "modal") return(modal)

  # The same labels the plotting methods use, so a posterior column and a plot
  # legend name the same class the same way.
  colnames(resp) <- paste("Class", seq_len(ncol(resp)))
  if (type == "posterior") return(resp)

  data.frame(class       = modal,
             probability = resp[cbind(seq_len(nrow(resp)), modal)],
             resp,
             check.names = FALSE)
}

# Weighted average posterior probability by modal class. The weights matter
# whenever a case stands for more than one observation, which is the norm for
# the frequency-weighted pattern files these models are often fitted to.
.ave_pp <- function(resp, w, K) {
  modal  <- max.col(resp, ties.method = "first")
  ave_pp <- matrix(NA_real_, nrow = K, ncol = K)
  for (k in seq_len(K)) {
    idx <- modal == k
    if (any(idx))
      ave_pp[k, ] <- colSums(resp[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
  }
  ave_pp
}

# ==============================================================================
# Internal helpers for summary.mixture_model
# ==============================================================================

# Bind collected table rows into one clean data frame (NULL when empty).
.rbind_tidy <- function(rows) {
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

# Format p-values to publication conventions.
# Values below .001 are shown as "< .001"; NaN/NA appear as a dash.
.fmt_pval <- function(p) {
  if (is.na(p) || is.nan(p)) return("       -")
  if (p < 0.001)              return("  < .001")
  sprintf("   %5.3f", p)
}

# Per-covariate omnibus Wald tests, printed under the class-predictor table.
#
# The coefficient table above answers "does this covariate distinguish class c
# from the reference class?" once per contrast. Two things are missing from it,
# and this table supplies both:
#
#   * Multiplicity. A 3-class model with five covariates prints twelve
#     coefficient tests, and some of them will be significant.
#   * For a covariate with more than two categories the table contains no test
#     of the covariate at all — a three-level marital status becomes two dummies
#     times two class contrasts, and none of those four rows answers "does
#     marital status predict class membership?", which is the question that was
#     asked. The omnibus test is the one that does.
#
# Nothing is printed when every covariate is a single column and there are only
# two classes, since each omnibus test is then the square of a z already shown.
.print_covariate_omnibus <- function(object, sm_sub, ref_class) {
  params <- sm_sub$parameters
  if (is.null(params$beta) || ncol(params$beta) == 0L) return(invisible(NULL))
  K <- object$n_components
  if (K < 2L) return(invisible(NULL))

  terms <- params$terms
  if (is.null(terms) || length(terms) != ncol(params$beta))
    terms <- colnames(params$beta) %||% paste0("V", seq_len(ncol(params$beta)))

  cov_terms <- unique(setdiff(terms, "Intercept"))
  if (!length(cov_terms)) return(invisible(NULL))
  if (K == 2L && length(cov_terms) == length(terms) - sum(terms == "Intercept"))
    return(invisible(NULL))

  rows <- lapply(cov_terms, function(tm) {
    cols <- which(terms == tm)
    st   <- tryCatch(.wald_omnibus_core(params, K, ref_class, cols),
                     error = function(e) NULL)
    if (is.null(st)) return(NULL)
    data.frame(term = tm, chi2 = st$W, df = st$df, p = st$p_value,
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(rows) || !nrow(rows)) return(invisible(NULL))

  cat("\nOMNIBUS TEST PER COVARIATE (effect across all classes)\n")
  cat("---------------------------------------------------------\n")
  # Same shortening policy as the coefficient table above, so a long variable
  # name reads identically in both places.
  disp    <- .shorten_labels(rows$term)
  label_w <- .label_width(disp, min = 20L)
  cat(sprintf("  %-*s %11s %4s  %s\n", label_w, "", "Wald Chi2", "df", "P-Value"))
  for (i in seq_len(nrow(rows)))
    cat(sprintf("  %-*s %11.3f %4d  %s\n",
                label_w, disp[i], rows$chi2[i], rows$df[i], .fmt_pval(rows$p[i])))
  .cat_label_legend(disp)
  # The Wald statistic is non-monotone in the effect size: a covariate strong
  # enough to nearly separate a class drives it back towards zero
  # (Hauck & Donner, 1977). The test gates in one direction only.
  cat("  Note: a non-significant test beside large coefficients can be the\n")
  cat("        Hauck-Donner effect; confirm with wald_omnibus_test().\n")
  invisible(rows)
}

# Omnibus Wald test for equality of K means (continuous distal outcomes).
#
# When Sigma_mu (the full K x K sandwich variance-covariance matrix of the
# means) is available, uses the full-covariance formulation:
#   W = c^T V^{-1} c,  where c = R * mu,  V = R * Sigma_mu * R^T
#   R = contrast matrix [class k vs class 1, k = 2..K]  (df = K-1)
#
# This is the robust Wald formulation of Bakk, Oberski and Vermunt
# (2014), accounting for the cross-class covariance induced by BCH weights.
#
# Falls back to the diagonal (precision-weighted) approximation when
# Sigma_mu is not stored (e.g. for non-BCH structural models).
#
# Returns a list with stat, df, and p.
.wald_omnibus_means <- function(means, ses, K, Sigma_mu = NULL) {
  if (K <= 1L) return(list(stat = NA_real_, df = 0L, p = NA_real_))
  df <- K - 1L

  if (!is.null(Sigma_mu) && all(is.finite(Sigma_mu))) {
    # Full sandwich Wald: contrast matrix R (K-1 x K), class 2..K vs class 1
    R     <- cbind(-1, diag(K - 1L))
    V     <- R %*% Sigma_mu %*% t(R)
    theta <- R %*% means
    W     <- tryCatch(
      as.numeric(t(theta) %*% solve(V) %*% theta),
      error = function(e) NA_real_
    )
  } else {
    # Fallback: diagonal precision-weighted approximation
    prec   <- 1 / pmax(ses^2, 1e-15)
    mu_bar <- sum(prec * means) / sum(prec)
    W      <- sum(prec * (means - mu_bar)^2)
  }

  p <- if (is.na(W)) NA_real_ else pchisq(W, df = df, lower.tail = FALSE)
  list(stat = W, df = df, p = p)
}

# Omnibus Wald test for equality of class effects in distal_pooled.
# Tests H0: beta[m, k] = beta[m, ref] for all k != ref and all m.
# df = (M - 1) * (K - 1).
.wald_omnibus_pooled <- function(beta_mat, Hessian, K, D_cov, ref_class) {
  M_minus_1 <- nrow(beta_mat)
  L         <- K + D_cov
  if (M_minus_1 == 0L || K <= 1L)
    return(list(stat = NA_real_, df = 0L, p = NA_real_))

  Sigma   <- pinv(-Hessian)
  non_ref <- setdiff(seq_len(K), ref_class)
  n_ctr   <- M_minus_1 * (K - 1L)
  R       <- matrix(0, nrow = n_ctr, ncol = M_minus_1 * L)

  row_i <- 1L
  for (m in seq_len(M_minus_1)) {
    for (k in non_ref) {
      R[row_i, (m - 1L) * L + k]         <-  1
      R[row_i, (m - 1L) * L + ref_class] <- -1
      row_i <- row_i + 1L
    }
  }

  # beta_mat is (M-1) x L; vectorise row-major to match Hessian block ordering
  beta_vec <- as.vector(t(beta_mat))
  r_vec    <- R %*% beta_vec
  V        <- R %*% Sigma %*% t(R)
  W        <- tryCatch(as.numeric(t(r_vec) %*% pinv(V) %*% r_vec),
                       error = function(e) NA_real_)
  p        <- if (is.na(W)) NA_real_ else pchisq(W, df = n_ctr, lower.tail = FALSE)
  list(stat = W, df = n_ctr, p = p)
}

# Predicted outcome probabilities for one class in a distal_pooled model,
# evaluated at covariates = 0 (i.e., the class intercept only).
.pred_probs_pooled <- function(beta_mat, K, D_cov, k) {
  L        <- K + D_cov
  U_k      <- matrix(0, nrow = 1, ncol = L)
  U_k[1, k] <- 1                                   # one-hot class indicator
  as.vector(distal_forward(U_k, beta_mat))
}

# Predicted outcome probabilities for one class in a distal_regression model,
# evaluated at covariates = 0 (i.e., using the intercept column only).
.pred_probs_reg <- function(betas_k, D) {
  U       <- matrix(0, nrow = 1, ncol = D)
  U[1, 1] <- 1                                     # intercept; covariates at 0
  as.vector(distal_forward(U, betas_k))
}

#' Summarise a Fitted Mixture Model
#'
#' @description
#' Prints a detailed summary of the structural model parameters. Depending on
#' which structural model was fitted, this includes covariate regression
#' coefficients (as odds ratios with 95\% confidence intervals and p-values),
#' distal outcome means, or class-specific regression effects. If no structural
#' model is present, a notice is printed directing the user to
#' \code{\link{measurement_summary}}.
#'
#' @param object A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ref_class Integer. The reference latent class for pairwise contrasts.
#'   Defaults to the first class (\code{1}).
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly, a list holding the printed numbers in tidy form, ready
#'   for further use: \code{$coefficients} (one row per class contrast and
#'   covariate: estimate, SE, p, odds ratio and its confidence limits),
#'   \code{$omnibus} (the per-covariate omnibus Wald tests), and, when a
#'   distal outcome is present, \code{$outcome} (predicted probabilities or
#'   class means/estimates with their tests). Returns \code{NULL} when the
#'   model has no structural part.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' summary(fit)
#'
#' @export
summary.mixture_model <- function(object, ref_class = NULL, ...) {
  K <- object$n_components
  if (is.null(ref_class)) ref_class <- 1

  # Input validation: ref_class must be a valid class index.
  # Without this guard the function starts printing output, then crashes
  # mid-way with a cryptic "subscript out of bounds" error.
  if (!is.numeric(ref_class) || length(ref_class) != 1 ||
      ref_class < 1 || ref_class > K)
    stop(sprintf(
      "ref_class must be an integer between 1 and %d. Got: %s",
      K, ref_class
    ))

  # Repeated ahead of everything else: a structural result read off a fit whose
  # measurement model has a collapsed class is not worth interpreting, so the
  # note has to come before the coefficients rather than after them.
  .print_degenerate_note(object)

  if (is.null(object$sm)) {
    cat("Notice: No structural model found. Use measurement_summary() for item parameters.\n")
    return(invisible())
  }

  # Everything printed below is also collected here and returned invisibly,
  # so vignettes and downstream code can use the numbers without re-deriving
  # them from the model internals.
  out <- list(ref_class  = ref_class,
              n_classes  = K,
              n_steps    = object$n_steps,
              correction = object$correction,
              assignment = object$assignment)

  cat("=========================================================\n")
  cat("             STRUCTURAL MODEL SUMMARY                    \n")
  cat("=========================================================\n")

  # A. Covariate model
  sm_sub <- NULL
  if (inherits(object$sm, "covariate")) sm_sub <- object$sm
  if (inherits(object$sm, "nested") && "predictor" %in% names(object$sm$models))
    sm_sub <- object$sm$models$predictor

  if (!is.null(sm_sub) && !is.null(sm_sub$parameters$hessian)) {
    cat("\nCATEGORICAL LATENT VARIABLE REGRESSION (CLASS PREDICTORS)\n")
    cat(sprintf("Reference Class: %d\n", ref_class))
    # A covariate standard error can come from any of several estimators that
    # differ by up to a factor of two, so the printed table names the one used
    # rather than leaving the reader to guess (Bakk, Oberski & Vermunt, 2014).
    cat(sprintf("Standard errors: %s\n",
                sm_sub$parameters$V_method %||%
                  if (!is.null(sm_sub$parameters$V_robust))
                    "Survey-robust (linearization)" else "Q-function Hessian"))
    cat("---------------------------------------------------------\n")

    betas     <- sm_sub$parameters$beta
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else paste0("V", 1:D)
    # The label column sizes itself to the longest (shortened) predictor name:
    # a dummy-coded predictor name runs past fifteen characters routinely, and
    # names past the shortening cap are abbreviated with a key printed under
    # the table. The data frame returned by summary() keeps the full names.
    disp      <- .shorten_labels(var_names)
    label_w   <- .label_width(disp, min = 20L)
    Sigma     <- if (!is.null(sm_sub$parameters$V_robust))
      sm_sub$parameters$V_robust else pinv(-sm_sub$parameters$hessian)

    cat(sprintf("  %-*s %9s  %s  %s\n", label_w,
                "", "OR", "       [95% CI]       ", "P-Value"))

    cov_rows <- vector("list", (K - 1L) * D)
    row_i    <- 0L
    for (c in setdiff(1:K, ref_class)) {
      cat(sprintf("\nClass %d ON\n", c))
      for (v in 1:D) {
        est     <- betas[c, v] - betas[ref_class, v]
        idx_c   <- (c - 1) * D + v
        idx_ref <- (ref_class - 1) * D + v
        var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
          2 * Sigma[idx_c, idx_ref]
        se    <- sqrt(max(0, var_diff))
        z_val <- est / se
        p_val <- 2 * (1 - pnorm(abs(z_val)))
        cat(sprintf("  %-*s %9.3f  [%9.3f, %9.3f]  %s\n",
                    label_w, disp[v], exp(est),
                    exp(est - 1.96 * se), exp(est + 1.96 * se),
                    .fmt_pval(p_val)))
        row_i <- row_i + 1L
        cov_rows[[row_i]] <- data.frame(
          class = c, term = var_names[v], estimate = est, se = se,
          z = z_val, p = p_val, OR = exp(est),
          OR_lower = exp(est - 1.96 * se), OR_upper = exp(est + 1.96 * se),
          stringsAsFactors = FALSE)
      }
    }
    out$coefficients <- .rbind_tidy(cov_rows[seq_len(row_i)])
    .cat_label_legend(disp)

    out$omnibus <- .print_covariate_omnibus(object, sm_sub, ref_class)
  }

  # B0. Categorical distal outcome with no covariate (distal_categorical).
  #     This section is checked before B because distal_categorical inherits
  #     from distal_pooled; without the explicit class check below, section B
  #     would match first and display the wrong header.
  cat_sub <- NULL
  if (inherits(object$sm, "distal_categorical"))
    cat_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_categorical"))
    cat_sub <- object$sm$models$distal

  if (!is.null(cat_sub) && !is.null(cat_sub$parameters$beta_pooled)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (CLASS PROBABILITIES)\n")
    cat("---------------------------------------------------------\n")

    beta_mat  <- cat_sub$parameters$beta_pooled
    M_minus_1 <- nrow(beta_mat)
    K_distal  <- K
    # distal_categorical has no covariate columns; D_cov is always 0.
    D_cov <- ncol(beta_mat) - K_distal
    M     <- M_minus_1 + 1L

    if (M_minus_1 == 0L) {
      cat("  (Constant outcome - no parameters to display)\n")
    } else {
      Sigma <- pinv(-cat_sub$parameters$hessian)

      # --- Omnibus test ---
      omni <- .wald_omnibus_pooled(beta_mat, cat_sub$parameters$hessian,
                                   K_distal, D_cov, ref_class)
      if (!is.na(omni$stat)) {
        cat(sprintf(
          "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
          omni$df, omni$stat, .fmt_pval(omni$p)))
      }

      # --- Predicted probabilities ---
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cat("\nPredicted Probabilities:\n")
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov = 0L, k)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table ---
      or_rows <- list()
      cat(sprintf("\nPairwise Odds Ratios (Reference: Class %d)\n", ref_class))
      cat("                     OR       [95% CI]        P-Value\n")
      for (m in seq_len(M_minus_1)) {
        cat(sprintf("\nOutcome Category %d (vs Category 1) ON\n", m + 1L))
        cat("  Latent Class:\n")
        for (c in setdiff(seq_len(K_distal), ref_class)) {
          est      <- beta_mat[m, c] - beta_mat[m, ref_class]
          idx_c    <- (m - 1L) * K_distal + c
          idx_ref  <- (m - 1L) * K_distal + ref_class
          var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
            2 * Sigma[idx_c, idx_ref]
          se    <- sqrt(max(0, var_diff))
          z_val <- if (se > 0) est / se else NA_real_
          p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
          cat(sprintf("    Class %d        %7.3f  [%6.3f, %6.3f]  %s\n",
                      c, exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se),
                      .fmt_pval(p_val)))
          or_rows[[length(or_rows) + 1L]] <- data.frame(
            category = m + 1L, class = c, estimate = est, se = se, p = p_val,
            OR = exp(est), OR_lower = exp(est - 1.96 * se),
            OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
        }
      }
      out$outcome <- list(
        type    = "categorical",
        omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
        predicted_probabilities = pp,
        odds_ratios = .rbind_tidy(or_rows))
    }
  }

  # B. Categorical distal outcome with a pooled covariate slope (distal_pooled).
  #    Excludes distal_categorical, which is handled in section B0 above.
  pooled_sub <- NULL
  if (inherits(object$sm, "distal_pooled") &&
      !inherits(object$sm, "distal_categorical")) pooled_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_pooled") &&
      !inherits(object$sm$models$distal, "distal_categorical"))
    pooled_sub <- object$sm$models$distal

  if (!is.null(pooled_sub) && !is.null(pooled_sub$parameters$beta_pooled)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (POOLED SLOPES)\n")
    cat("---------------------------------------------------------\n")

    beta_mat  <- pooled_sub$parameters$beta_pooled
    M_minus_1 <- nrow(beta_mat)
    K_distal  <- K
    D_cov     <- ncol(beta_mat) - K_distal
    M         <- M_minus_1 + 1L
    Sigma     <- pinv(-pooled_sub$parameters$hessian)
    var_names <- if (D_cov > 0) paste0("Z", seq_len(D_cov)) else character(0)

    if (M_minus_1 > 0L) {

      # --- Omnibus test ---
      omni <- .wald_omnibus_pooled(beta_mat, pooled_sub$parameters$hessian,
                                   K_distal, D_cov, ref_class)
      if (!is.na(omni$stat)) {
        cat(sprintf(
          "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
          omni$df, omni$stat, .fmt_pval(omni$p)))
      }

      # --- Predicted probabilities (primary display) ---
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cov_note <- if (D_cov > 0) " (covariates held at zero)" else ""
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        probs <- .pred_probs_pooled(beta_mat, K_distal, D_cov, k)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Pairwise OR table (secondary) ---
      or_rows  <- list()
      cov_rows <- list()
      cat(sprintf("\nPairwise Odds Ratios (Reference: Class %d)\n", ref_class))
      cat("                     OR       [95% CI]        P-Value\n")
      for (m in seq_len(M_minus_1)) {
        cat(sprintf("\nOutcome Category %d (vs Category 1) ON\n", m + 1L))
        cat("  Latent Class:\n")
        for (c in setdiff(seq_len(K_distal), ref_class)) {
          est      <- beta_mat[m, c] - beta_mat[m, ref_class]
          idx_c    <- (m - 1L) * (K_distal + D_cov) + c
          idx_ref  <- (m - 1L) * (K_distal + D_cov) + ref_class
          var_diff <- Sigma[idx_c, idx_c] + Sigma[idx_ref, idx_ref] -
            2 * Sigma[idx_c, idx_ref]
          se    <- sqrt(max(0, var_diff))
          z_val <- if (se > 0) est / se else NA_real_
          p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
          cat(sprintf("    Class %d        %7.3f  [%6.3f, %6.3f]  %s\n",
                      c, exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se),
                      .fmt_pval(p_val)))
          or_rows[[length(or_rows) + 1L]] <- data.frame(
            category = m + 1L, class = c, estimate = est, se = se, p = p_val,
            OR = exp(est), OR_lower = exp(est - 1.96 * se),
            OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
        }
        if (D_cov > 0) {
          cat("  Covariates (Pooled Slope):\n")
          for (v in seq_len(D_cov)) {
            est   <- beta_mat[m, K_distal + v]
            idx   <- (m - 1L) * (K_distal + D_cov) + K_distal + v
            se    <- sqrt(max(0, Sigma[idx, idx]))
            z_val <- if (se > 0) est / se else NA_real_
            p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
            cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]  %s\n",
                        var_names[v], exp(est),
                        exp(est - 1.96 * se), exp(est + 1.96 * se),
                        .fmt_pval(p_val)))
            cov_rows[[length(cov_rows) + 1L]] <- data.frame(
              category = m + 1L, term = var_names[v], estimate = est, se = se,
              p = p_val, OR = exp(est), OR_lower = exp(est - 1.96 * se),
              OR_upper = exp(est + 1.96 * se), stringsAsFactors = FALSE)
          }
        }
      }
      out$outcome <- list(
        type    = "categorical",
        omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
        predicted_probabilities = pp,
        odds_ratios = .rbind_tidy(or_rows),
        covariate_effects = .rbind_tidy(cov_rows))
    }
  }

  # C. Distal regression (moderated)
  distal_sub <- NULL
  if (inherits(object$sm, "distal_regression")) distal_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_regression"))
    distal_sub <- object$sm$models$distal

  if (!is.null(distal_sub) && !is.null(distal_sub$parameters$betas)) {
    cat("\nCATEGORICAL DISTAL OUTCOME (CLASS-SPECIFIC SLOPES)\n")
    cat("---------------------------------------------------------\n")

    distal_betas <- distal_sub$parameters$betas
    K_distal     <- dim(distal_betas)[1]
    M_minus_1    <- dim(distal_betas)[2]
    D_distal     <- dim(distal_betas)[3]
    M            <- M_minus_1 + 1L

    if (M_minus_1 == 0L) {
      cat("  (Constant outcome - no parameters to display)\n")
    } else {
      var_names <- c("Intercept", paste0("Z", seq_len(D_distal - 1L)))
      cov_note  <- if (D_distal > 1L) " (covariates held at zero)" else ""

      # --- Predicted probabilities (primary display) ---
      pp <- matrix(NA_real_, K_distal, M,
                   dimnames = list(paste("Class", seq_len(K_distal)),
                                   paste("Cat", seq_len(M))))
      cat(sprintf("\nPredicted Probabilities%s:\n", cov_note))
      cat(sprintf("  %-12s", ""))
      for (m in seq_len(M)) cat(sprintf(" Cat %-4d", m))
      cat("\n")
      for (k in seq_len(K_distal)) {
        betas_k <- matrix(distal_betas[k, , ], nrow = M_minus_1, ncol = D_distal)
        probs   <- .pred_probs_reg(betas_k, D_distal)
        pp[k, ] <- probs
        cat(sprintf("  Class %-6d", k))
        for (m in seq_len(M)) cat(sprintf("  %6.3f ", probs[m]))
        cat("\n")
      }

      # --- Class-specific OR tables (secondary) ---
      est_rows <- list()
      cat("\nClass-Specific Estimates\n")
      cat("                     OR       [95% CI]        P-Value\n")
      for (k in seq_len(K_distal)) {
        cat(sprintf("\nClass %d:\n", k))
        Sigma <- if (!is.null(distal_sub$parameters$hessians) &&
                     length(distal_sub$parameters$hessians) >= k)
          pinv(-distal_sub$parameters$hessians[[k]])
        else
          matrix(0, M_minus_1 * D_distal, M_minus_1 * D_distal)

        for (m in seq_len(M_minus_1)) {
          cat(sprintf("  Outcome Category %d (vs Category 1) ON\n", m + 1L))
          for (v in seq_len(D_distal)) {
            est   <- distal_betas[k, m, v]
            idx   <- (m - 1L) * D_distal + v
            se    <- sqrt(max(0, Sigma[idx, idx]))
            z_val <- if (se > 0) est / se else NA_real_
            p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
            if (se > 0) {
              cat(sprintf("    %-13s %7.3f  [%6.3f, %6.3f]  %s\n",
                          var_names[v], exp(est),
                          exp(est - 1.96 * se), exp(est + 1.96 * se),
                          .fmt_pval(p_val)))
            } else {
              cat(sprintf("    %-13s %7.3f  [   N/A,    N/A]       N/A\n",
                          var_names[v], exp(est)))
            }
            est_rows[[length(est_rows) + 1L]] <- data.frame(
              class = k, category = m + 1L, term = var_names[v],
              estimate = est, se = if (se > 0) se else NA_real_, p = p_val,
              OR = exp(est), stringsAsFactors = FALSE)
          }
        }
      }
      out$outcome <- list(
        type = "categorical",
        predicted_probabilities = pp,
        estimates = .rbind_tidy(est_rows))
    }
  }

  # D. Continuous distal (means)
  cont_sub <- NULL
  if (inherits(object$sm, "distal_continuous")) cont_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous"))
    cont_sub <- object$sm$models$distal

  if (!is.null(cont_sub) && !is.null(cont_sub$parameters$means)) {
    cat("\nCONTINUOUS DISTAL OUTCOME (MEANS)\n")
    cat("---------------------------------------------------------\n")

    means    <- as.vector(cont_sub$parameters$means)
    ses      <- as.vector(cont_sub$parameters$ses)
    Sigma_mu <- cont_sub$parameters$Sigma_mu   # NULL for non-BCH models

    # Omnibus Wald test for equality of class means
    omni <- .wald_omnibus_means(means, ses, K, Sigma_mu = Sigma_mu)
    if (!is.na(omni$stat)) {
      cat(sprintf(
        "\nOmnibus test (class differences): Wald chi^2(%d) = %.2f, p%s\n",
        omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    cat("\n                 Mean       [95% CI]        SE\n")
    for (k in seq_len(K)) {
      mu <- means[k]
      se <- ses[k]
      cat(sprintf("  Class %d      %7.3f  [%6.3f, %6.3f]   %7.3f\n",
                  k, mu, mu - 1.96 * se, mu + 1.96 * se, se))
    }
    # Which classes differ, not just whether any do. The omnibus above is one
    # joint statistic; what gets written up is the contrast, and computing it
    # by hand from the two standard errors gets it wrong, because the class
    # means are correlated. The categorical outcome sections print their own
    # contrasts already (the odds-ratio tables), so this is the branch that had
    # nothing.
    #
    # Printed only for a handful of classes: at K = 8 the all-pairs table is
    # twenty-eight rows and would bury the summary it is part of. Above the cap
    # the table is a call away, and either way it is returned.
    contrasts <- tryCatch(outcome_contrasts(object), error = function(e) NULL)
    if (!is.null(contrasts)) {
      if (K <= 5L) {
        cat("\nPairwise class differences:\n")
        cat("                    Difference       [95% CI]        P-Value\n")
        for (i in seq_len(nrow(contrasts))) {
          cat(sprintf(
            "  Class %d vs %d       %7.3f  [%6.3f, %6.3f]  %s\n",
            contrasts$class[i], contrasts$reference[i], contrasts$estimate[i],
            contrasts$lower[i], contrasts$upper[i], .fmt_pval(contrasts$p[i])))
        }
      } else {
        cat(sprintf(paste0("\n%d pairwise class differences are not printed at ",
                           "%d classes; see outcome_contrasts(fit).\n"),
                    nrow(contrasts), K))
      }
    }

    out$outcome <- list(
      type    = "continuous",
      omnibus = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      means   = data.frame(class = seq_len(K), mean = means, se = ses,
                           lower = means - 1.96 * ses,
                           upper = means + 1.96 * ses),
      contrasts = contrasts)
  }

  # E. Continuous distal regression
  cont_reg_sub <- NULL
  if (inherits(object$sm, "distal_continuous_regression")) cont_reg_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_regression"))
    cont_reg_sub <- object$sm$models$distal

  if (!is.null(cont_reg_sub) && !is.null(cont_reg_sub$parameters$betas)) {
    cat("\nCONTINUOUS DISTAL REGRESSION (Y ~ Z * Class)\n")
    cat("---------------------------------------------------------\n")

    betas     <- cont_reg_sub$parameters$betas
    ses       <- cont_reg_sub$parameters$ses
    D         <- ncol(betas)
    var_names <- if (!is.null(colnames(betas))) colnames(betas) else
      c("Intercept", paste0("Z", seq_len(D - 1L)))
    disp      <- .shorten_labels(var_names)
    label_w   <- .label_width(disp, min = 13L)
    hdr_pad   <- strrep(" ", 17L + label_w - 13L)

    # Omnibus Wald test on class intercepts (class-specific means at Z = 0).
    # Uses the model-based SE for each intercept (sigma^2 * B_inv_k[1,1]),
    # assuming classes are independent (a Wald test of equality).
    intercepts <- betas[, 1L]
    int_ses    <- ses[, 1L]   # already model-based if fitted with BCH v2
    prec_int   <- 1 / pmax(int_ses^2, 1e-15)
    mu_bar_int <- sum(prec_int * intercepts) / sum(prec_int)
    W_stat_int <- sum(prec_int * (intercepts - mu_bar_int)^2)
    omni <- list(stat = W_stat_int, df = K - 1L,
                 p = pchisq(W_stat_int, df = K - 1L, lower.tail = FALSE))
    if (!is.na(omni$stat)) {
      cov_note <- if (D > 1L) " (at covariate zero)" else ""
      cat(sprintf(
        "\nOmnibus test (class differences%s): Wald chi^2(%d) = %.2f, p%s\n",
        cov_note, omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    est_rows <- list()
    cat("\n")
    for (k in seq_len(K)) {
      cat(sprintf("Class %d:\n", k))
      cat(hdr_pad, "Estimate   [95% CI]        P-Value\n", sep = "")
      for (v in seq_len(D)) {
        est   <- betas[k, v]
        se    <- ses[k, v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("  %-*s %7.3f  [%6.3f, %6.3f]  %s\n",
                    label_w, disp[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
        est_rows[[length(est_rows) + 1L]] <- data.frame(
          class = k, term = var_names[v], estimate = est, se = se, p = p_val,
          lower = est - 1.96 * se, upper = est + 1.96 * se,
          stringsAsFactors = FALSE)
      }
      cat("\n")
    }

    #  Per-covariate Wald(=) tests: H0: slope_k equal across all classes
    # Uses the diagonal independence approximation (separate per-class
    # regressions).
    # Contrast matrix R = [-1 | I_{K-1}], df = K-1.
    eq_rows <- list()
    if (K > 1L && D > 1L) {
      cat("---------------------------------------------------------\n")
      cat("Wald tests (equality of slopes across classes):\n")
      cat(sprintf("  %-*s   Wald(%s)%s  P-Value\n",
                  label_w, "", paste0("chi^2(", K - 1L, ")"), ""))
      R_eq <- cbind(-1, diag(K - 1L))
      for (v in 2L:D) {
        theta_v <- betas[, v]
        var_v   <- ses[, v]^2
        V_c     <- R_eq %*% diag(var_v) %*% t(R_eq)
        th_c    <- R_eq %*% theta_v
        W_v     <- tryCatch(
          as.numeric(t(th_c) %*% solve(V_c) %*% th_c),
          error = function(e) NA_real_)
        p_v     <- if (!is.na(W_v))
          pchisq(W_v, df = K - 1L, lower.tail = FALSE) else NA_real_
        cat(sprintf("  %-*s   %8.2f          %s\n",
                    label_w, disp[v], W_v, .fmt_pval(p_v)))
        eq_rows[[length(eq_rows) + 1L]] <- data.frame(
          term = var_names[v], chi2 = W_v, df = K - 1L, p = p_v,
          stringsAsFactors = FALSE)
      }
    }
    .cat_label_legend(disp)
    out$outcome <- list(
      type      = "continuous",
      omnibus   = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      estimates = .rbind_tidy(est_rows),
      slope_equality = .rbind_tidy(eq_rows))
  }

  # F. Continuous distal pooled
  cont_pool_sub <- NULL
  if (inherits(object$sm, "distal_continuous_pooled")) cont_pool_sub <- object$sm
  if (inherits(object$sm, "nested") && "distal" %in% names(object$sm$models) &&
      inherits(object$sm$models$distal, "distal_continuous_pooled"))
    cont_pool_sub <- object$sm$models$distal

  if (!is.null(cont_pool_sub) && !is.null(cont_pool_sub$parameters$beta_pooled)) {
    cat("\nCONTINUOUS DISTAL POOLED REGRESSION (Main Effects)\n")
    cat("---------------------------------------------------------\n")

    theta     <- as.vector(cont_pool_sub$parameters$beta_pooled)
    ses       <- as.vector(cont_pool_sub$parameters$ses)
    K_distal  <- K
    D_cov     <- length(theta) - K_distal
    # Use stored column names when available; fall back to Z1, Z2, ...
    stored_names <- colnames(cont_pool_sub$parameters$beta_pooled)
    var_names <- if (!is.null(stored_names) && length(stored_names) > K_distal)
      stored_names[(K_distal + 1L):length(stored_names)]
    else if (D_cov > 0) paste0("Z", seq_len(D_cov))
    else character(0)

    # Omnibus Wald test on class intercepts
    # When cov_theta is available (BCH step stored it), use the full
    # model-based contrast Wald: H0: int_k = int_1 for all k != 1.
    # This omnibus formulation accounts for the
    # covariance between intercept estimates.
    cov_theta  <- cont_pool_sub$parameters$cov_theta
    intercepts <- theta[seq_len(K_distal)]
    int_ses    <- ses[seq_len(K_distal)]

    if (!is.null(cov_theta) && all(is.finite(cov_theta))) {
      cov_int   <- cov_theta[seq_len(K_distal), seq_len(K_distal)]
      R_int     <- cbind(-1, diag(K_distal - 1L))
      V_contr   <- R_int %*% cov_int %*% t(R_int)
      theta_c   <- R_int %*% intercepts
      W_stat    <- tryCatch(
        as.numeric(t(theta_c) %*% solve(V_contr) %*% theta_c),
        error = function(e) NA_real_
      )
      omni <- list(stat = W_stat, df = K_distal - 1L,
                   p = if (is.na(W_stat)) NA_real_
                   else pchisq(W_stat, df = K_distal - 1L, lower.tail = FALSE))
    } else {
      omni <- .wald_omnibus_means(intercepts, int_ses, K_distal)
    }

    if (!is.na(omni$stat)) {
      cov_note <- if (D_cov > 0) " (at covariate zero)" else ""
      cat(sprintf(
        "\nOmnibus test (class differences%s): Wald chi^2(%d) = %.2f, p%s\n",
        cov_note, omni$df, omni$stat, .fmt_pval(omni$p)))
    }

    int_rows <- list()
    cat("\n  Latent Class (Intercepts):\n")
    cat("                 Estimate   [95% CI]        P-Value\n")
    for (k in seq_len(K_distal)) {
      est   <- theta[k]
      se    <- ses[k]
      z_val <- if (se > 0) est / se else NA_real_
      p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
      cat(sprintf("    Class %d      %7.3f  [%6.3f, %6.3f]  %s\n",
                  k, est, est - 1.96 * se, est + 1.96 * se,
                  .fmt_pval(p_val)))
      int_rows[[length(int_rows) + 1L]] <- data.frame(
        class = k, estimate = est, se = se, p = p_val,
        lower = est - 1.96 * se, upper = est + 1.96 * se,
        stringsAsFactors = FALSE)
    }

    slope_rows <- list()
    if (D_cov > 0) {
      disp    <- .shorten_labels(var_names)
      label_w <- .label_width(disp, min = 11L)
      cat("\n  Covariates (Pooled Slopes):\n")
      cat(strrep(" ", 17L + label_w - 11L),
          "Estimate   [95% CI]        P-Value\n", sep = "")
      for (v in seq_len(D_cov)) {
        est   <- theta[K_distal + v]
        se    <- ses[K_distal + v]
        z_val <- if (se > 0) est / se else NA_real_
        p_val <- if (!is.na(z_val)) 2 * (1 - pnorm(abs(z_val))) else NA_real_
        cat(sprintf("    %-*s %7.3f  [%6.3f, %6.3f]  %s\n",
                    label_w, disp[v], est, est - 1.96 * se, est + 1.96 * se,
                    .fmt_pval(p_val)))
        slope_rows[[length(slope_rows) + 1L]] <- data.frame(
          term = var_names[v], estimate = est, se = se, p = p_val,
          lower = est - 1.96 * se, upper = est + 1.96 * se,
          stringsAsFactors = FALSE)
      }
      .cat_label_legend(disp)
    }
    out$outcome <- list(
      type       = "continuous",
      omnibus    = data.frame(chi2 = omni$stat, df = omni$df, p = omni$p),
      intercepts = .rbind_tidy(int_rows),
      covariate_effects = .rbind_tidy(slope_rows))
  }

  cat("=========================================================\n")
  invisible(out)
}

#' Fit a Latent Mixture Model (LCA / LPA)
#'
#' @description
#' The core estimation function. Fits a latent class analysis (LCA) or latent
#' profile analysis (LPA) model using the EM algorithm. Optionally fits a
#' structural model (covariates or distal outcomes) using 1-, 2-, or 3-step
#' estimation with optional bias correction.
#'
#' @param X A numeric matrix or data frame of indicator variables for the
#'   measurement model. Rows are observations; columns are items or variables.
#' @param Y Optional numeric matrix or data frame of outcome or covariate
#'   variables for the structural model. Must be provided when
#'   \code{structural} is not \code{NULL}.
#' @param n_components Positive integer. Number of latent classes (or profiles)
#'   to estimate. Default is \code{2}.
#' @param measurement Character string or named list specifying the measurement
#'   model type. Accepted strings: \code{"binary"} / \code{"bernoulli"},
#'   \code{"categorical"} / \code{"multinoulli"},
#'   \code{"continuous"} / \code{"gaussian_diag"},
#'   \code{"gaussian"} / \code{"gaussian_unit"},
#'   \code{"count"} / \code{"poisson"}.
#'   Missing values are handled automatically: any indicator column containing
#'   \code{NA} is estimated with a full-information (FIML) variant that masks the
#'   missing cells under a missing-at-random assumption, while complete columns
#'   use the faster complete-data estimator. A single specification (e.g.
#'   \code{"binary"}) therefore covers both complete and incomplete data, and the
#'   fitted object reports any missingness it found. Cases missing on
#'   \emph{every} indicator are the exception: they carry no information about
#'   class membership, so they are deleted before estimation and
#'   reported by \code{print()} and in \code{$missing_data$n_empty_rows}. The explicit \code{"*_nan"}
#'   forms (e.g. \code{"binary_nan"}, \code{"continuous_nan"}) remain accepted as
#'   aliases that force the missing-data variant.
#'   Pass a named list to specify a mixed (nested) measurement model with
#'   different variable types; each block's missing-data handling is resolved
#'   from the columns it governs. Default is \code{"binary"}.
#' @param structural Character string specifying the structural model type.
#'   One of \code{"covariate"}, \code{"distal_regression"},
#'   \code{"distal_pooled"}, \code{"distal_continuous"},
#'   \code{"distal_continuous_regression"}. Requires \code{Y}. Default is
#'   \code{NULL} (measurement model only).
#' @param n_steps Integer. Estimation approach: \code{1} for simultaneous
#'   1-step, \code{2} for 2-step, or \code{3} for bias-corrected 3-step.
#'   Default is \code{1}.
#' @param correction Character. Bias correction for 3-step estimation.
#'   One of \code{"none"}, \code{"BCH"}, or \code{"ML"}. Ignored when
#'   \code{n_steps} is not \code{3}. Default is \code{"none"}.
#' @param n_init Positive integer. Number of random restarts. The solution
#'   with the highest log-likelihood is retained. Default is \code{20}.
#' @param max_iter Positive integer. Maximum EM iterations per restart.
#'   Default is \code{1000}.
#' @param random_state Optional integer seed for reproducibility. Default is
#'   \code{NULL}.
#' @param order_by_size Logical. If \code{TRUE} (default), classes are sorted
#'   from largest to smallest after fitting.
#' @param weights Optional numeric vector of length \code{nrow(X)} for survey
#'   or case weights. Default is \code{NULL} (equal weights of 1).
#' @param weight_type Either \code{"sampling"} (default; weights are rescaled to
#'   sum to the number of cases) or \code{"frequency"} (weights are counts of
#'   identical cases and set the effective sample size).
#' @param strata Optional vector of stratum identifiers for complex survey designs.
#' @param cluster Optional vector of cluster identifiers for complex survey designs.
#' @param bayes_constants Optional named list of prior strengths
#'   (\code{latent}, \code{categorical}, \code{poisson}, \code{variances}), each
#'   defaulting to \code{1}. See \code{\link{fit_mixture}}.
#' @param refine Logical. If \code{TRUE} (default), applies L-BFGS refinement
#'   after EM convergence to optimize the penalized maximum likelihood.
#' @param warm_start Optional function of \code{(model_state, X, Y)} returning a
#'   starting model state for EM, or \code{NULL} to skip that start. Used by the
#'   group-varying measurement search to seed each fit from the pooled solution.
#'   \code{NULL} (default) uses only the usual random initializations.
#' @param se Character. How standard errors for a covariate (class-prediction)
#'   structural model are computed when \code{n_steps} is \code{2} or \code{3}.
#'   \code{"corrected"} (default) is the first-order corrected estimator of
#'   Bakk et al. (2014): the step-3 sandwich plus the variance
#'   propagated from step 1. \code{"robust"} keeps only the sandwich.
#'   \code{"hessian"} inverts the
#'   step-3 observed information alone. See \code{\link{covariate_se}} for the
#'   differences and when they matter. Ignored for other structural models and
#'   for \code{n_steps = 1}.
#' @param ... Additional arguments passed to the measurement or structural
#'   model constructors (e.g., \code{max_val} for multinoulli models).
#'
#' @return An object of class `mixture_model`, a list with:
#'   * `n_components` Number of latent classes.
#'   * `weights` Numeric vector of estimated class proportions.
#'   * `mm` Fitted measurement model state object.
#'   * `sm` Fitted structural model state object, or `NULL`.
#'   * `metrics` Named list: `ll` (log-likelihood), `aic`, `bic`, `sabic`,
#'     `n_params`, and `entropy` (relative entropy, 0-1 scale).
#'   * `log_resp` Matrix of log posterior class probabilities (n x K).
#'     Use `exp(fit$log_resp)` to obtain posterior probabilities.
#'   * `converged` Logical. Whether the EM algorithm converged.
#'   * `n_iter` Integer. Number of EM iterations run.
#'   * `step1_metrics` Named list of Step-1 fit indices (only when
#'     `n_steps = 3`).
#'
#' @examples
#' # Binary LCA with 3 classes and 5 random restarts
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 3, measurement = "binary", n_init = 5)
#' print(fit)
#' summary(fit)
#' measurement_summary(fit)
#'
#' # Continuous LPA (2 classes)
#' X_cont <- matrix(rnorm(300), nrow = 100)
#' fit_lpa <- fit_mixture(X_cont, n_components = 2, measurement = "continuous")
#'
#' # 3-step LCA with a covariate and ML correction
#' Z <- matrix(rnorm(100), nrow = 100)
#' fit_cov <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                        structural = "covariate",
#'                        n_steps = 3, correction = "ML", n_init = 5)
#' summary(fit_cov)
#'
#' @export
#' @importFrom stats complete.cases cov dnorm dpois optim pchisq plogis pnorm qlogis qnorm rbinom rnorm rpois runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
fit_mixture_internal <- function(X, Y = NULL, n_components = 2,
                                 measurement = "binary", structural = NULL,
                                 n_steps = 1, correction = "none",
                                 assignment = c("proportional", "modal"),
                                 n_init = 20,
                                 max_iter = 1000, random_state = NULL,
                                 order_by_size = TRUE, weights = NULL,
                                 weight_type = c("sampling", "frequency"),
                                 strata = NULL, cluster = NULL,
                                 refine = TRUE,
                                 bayes_constants = NULL,
                                 warm_start = NULL,
                                 se = c("corrected", "robust", "hessian"), ...) {

  weight_type <- match.arg(weight_type)
  se          <- match.arg(se)
  assignment  <- match.arg(assignment)
  bayes_constants <- .resolve_bayes_constants(bayes_constants)

  if (is.data.frame(X)) X <- as.matrix(X)
  # Convert Y through prepare_covariates() so that:
  #   - numeric columns are passed through unchanged
  #   - factor / character columns are dummy-coded (first level = reference)
  #   - column names are always preserved for display in summary()
  if (!is.null(Y)) Y <- prepare_covariates(Y)

  # --- Cases with no observed indicator data ----------------------------------
  # Dropped before anything else is computed, so that the sample size, the
  # missingness summary, the entropy, and the modal class counts all describe
  # the cases actually analysed. See .empty_rows() for why these cases distort
  # those quantities while contributing nothing to the likelihood.
  empty_rows <- .empty_rows(X)
  n_input_rows <- nrow(X)

  if (length(empty_rows) > 0L) {
    if (length(empty_rows) == n_input_rows)
      stop("Every case is missing on all indicators, so there are no data to ",
           "fit. This usually means the indicator columns were not read ",
           "correctly.", call. = FALSE)

    # Row-aligned arguments are length-checked against the input before any
    # subsetting: a mis-specified vector must still produce its own clear error
    # rather than being silently truncated to the retained rows.
    if (!is.null(weights) && length(weights) != n_input_rows)
      stop(sprintf("`weights` must have one entry per case (expected %d, got %d).",
                   n_input_rows, length(weights)), call. = FALSE)
    if (!is.null(strata) && length(strata) != n_input_rows)
      stop("Length of strata must match rows of X.", call. = FALSE)
    if (!is.null(cluster) && length(cluster) != n_input_rows)
      stop("Length of cluster must match rows of X.", call. = FALSE)

    keep <- setdiff(seq_len(n_input_rows), empty_rows)
    X    <- X[keep, , drop = FALSE]
    if (!is.null(Y))       Y       <- Y[keep, , drop = FALSE]
    if (!is.null(weights)) weights <- weights[keep]
    if (!is.null(strata))  strata  <- strata[keep]
    if (!is.null(cluster)) cluster <- cluster[keep]

    warning(sprintf(
      paste0("%d case%s had no observed value on any indicator and %s removed ",
             "before estimation (n = %d analysed). Rows: %s."),
      length(empty_rows), if (length(empty_rows) == 1L) "" else "s",
      if (length(empty_rows) == 1L) "was" else "were", length(keep),
      .abbreviate_indices(empty_rows)), call. = FALSE)
  }

  n_samples <- nrow(X)

  # --- Input Validation ---

  # n_components must be a positive integer. Values of 0 or below would either
  # produce degenerate output silently or crash with uninformative messages.
  if (!is.numeric(n_components) || length(n_components) != 1 ||
      !is.finite(n_components) || n_components < 1L)
    stop(sprintf(
      "n_components must be a positive integer (>= 1). Got: %s",
      n_components
    ))

  # n_steps must be 1, 2, or 3
  if (!n_steps %in% c(1L, 2L, 3L))
    stop(sprintf("n_steps must be 1, 2, or 3. Got: %d", n_steps))

  # correction must be one of the three supported values.
  valid_corrections <- c("none", "ML", "BCH")
  if (!correction %in% valid_corrections)
    stop(sprintf(
      "correction '%s' not recognized. Choose from: %s",
      correction, paste(valid_corrections, collapse = ", ")
    ))

  # Missing values in the indicator matrix are handled automatically. The
  # measurement descriptor is resolved against the data so complete columns keep
  # the fast complete-data estimator while columns containing NA switch to the
  # FIML variant that masks missing cells (MAR assumption). The user-facing
  # specification (e.g. "binary") therefore covers both complete and incomplete
  # data; explicit "*_nan" strings remain accepted as aliases.
  measurement_engine <- .resolve_emission_descriptor(measurement, X)

  # Summarise missingness so it can be reported in the fitted object's print and
  # measurement summaries, and so the estimator used is explicit downstream.
  item_missing <- colSums(is.na(X))

  # Structural-side missingness (covariates / distal outcomes in Y). Covariates
  # are completed under the class-invariant Gaussian marginal
  # (endogenous-constrained-x; Sterba, 2014); a missing distal outcome is
  # handled by FIML inside the structural likelihood.
  y_any_missing <- !is.null(Y) && anyNA(Y)
  y_n_missing   <- if (!is.null(Y)) sum(is.na(Y)) else 0L
  struct_handled <- if (y_any_missing) {
    if (!is.null(structural) && structural %in%
        c("covariate", "predict_class"))
      "endogenous-constrained-x (Sterba, 2014)"
    else
      "endogenous-constrained-x covariates; FIML outcome"
  } else NA_character_

  missing_data <- list(
    any_missing      = anyNA(X) || y_any_missing,
    n_missing        = sum(is.na(X)),
    n_cells          = length(X),
    prop_missing     = if (length(X) > 0) mean(is.na(X)) else 0,
    per_item         = item_missing,
    n_items_affected = sum(item_missing > 0),
    handled_by       = if (anyNA(X)) "FIML (MAR assumption)" else NA_character_,
    y_any_missing    = y_any_missing,
    y_n_missing      = y_n_missing,
    structural_handled_by = struct_handled,
    # Cases removed for having no observed indicator at all. Kept separate from
    # the FIML summary above, which describes only the analysed cases.
    n_empty_rows     = length(empty_rows),
    empty_rows       = empty_rows,
    n_input_rows     = n_input_rows
  )

  # Validate binary data when a Bernoulli family is requested.
  if (is.character(measurement) &&
      measurement %in% c("binary", "bernoulli", "binary_nan", "bernoulli_nan")) {
    # Two-valued indicators of any coding are converted on the way in (see
    # .recode_binary()), so anything reaching here has three or more distinct
    # values and is not dichotomous at all. The right move is a different
    # measurement model, which the message now says.
    valid_vals <- X[!is.na(X)]
    if (length(valid_vals) > 0 && !all(valid_vals %in% c(0, 1)))
      stop(sprintf(
        paste0("measurement = '%s' needs indicators with two values, and at ",
               "least one has more. Two-valued items are converted to 0/1 ",
               "automatically whatever their coding, so this is not a coding ",
               "problem: use measurement = \"categorical\" for ordered or ",
               "nominal items with three or more categories, or \"continuous\" ",
               "for a scale. Values found beyond {0, 1}: %s"),
        measurement,
        # head() rather than [1:5]: the old form indexed by the *count* of
        # offending cells, not the number of distinct offending values, so a
        # column with two bad values repeated many times printed "2, 3, NA, NA".
        paste(utils::head(sort(unique(valid_vals[!valid_vals %in% c(0, 1)])), 5L),
              collapse = ", ")
      ))
  }

  # Validate count data when a Poisson family is requested. Negative or
  # fractional values give dpois() a probability of zero for every class, so the
  # fit would fail with a log-likelihood of -Inf rather than a usable message.
  if (is.character(measurement) &&
      measurement %in% c("count", "poisson", "count_nan", "poisson_nan")) {
    valid_vals <- X[!is.na(X)]
    bad <- valid_vals[valid_vals < 0 | abs(valid_vals - round(valid_vals)) > 1e-8]
    if (length(bad) > 0)
      stop(sprintf(
        paste0("measurement = '%s' requires non-negative integer counts. ",
               "Found values outside this set: %s"),
        measurement,
        paste(utils::head(sort(unique(bad)), 5), collapse = ", ")
      ), call. = FALSE)
  }

  # Weights that are whole numbers adding up to far more than the number of rows
  # are almost always frequency counts, and treating them as sampling weights
  # would report a sample size of 31 patterns where the study had 631 people.
  # Say so rather than silently producing a BIC on the wrong scale.
  if (weight_type == "sampling" && .looks_like_frequencies(weights, n_samples))
    message("These weights look like frequency counts (whole numbers summing ",
            "to ", format(sum(weights)), " across ", n_samples, " rows). ",
            "They are being treated as sampling weights; use ",
            "weight_type = \"frequency\" if each row stands for that many cases.")

  wt <- .resolve_weights(weights, n_samples, weight_type)
  weights <- wt$weights
  n_eff   <- wt$n_eff

  # Survey design variables, when supplied, must align with the rows of X.
  # A design is considered present if either strata or cluster is given; the
  # other defaults so that every observation forms its own PSU or single
  # stratum, which leaves the linearization variance well defined.
  if (!is.null(strata) && length(strata) != n_samples)
    stop("Length of strata must match rows of X.")
  if (!is.null(cluster) && length(cluster) != n_samples)
    stop("Length of cluster must match rows of X.")
  has_survey_design <- !is.null(strata) || !is.null(cluster)

  # Structural model requires Y. Without this guard the SM is built but never
  # fitted (m_step_core gates on !is.null(Y)), so parameters$beta stays NULL
  # and every downstream function (coef, confint, wald tests) crashes with a
  # cryptic error rather than pointing here.
  if (!is.null(structural) && is.null(Y))
    stop(paste(
      "A structural model was specified but Y is NULL.",
      "Provide a Y matrix containing the outcome/covariate data,",
      "or set structural = NULL for a measurement-only model."
    ))

  model_state <- list(
    n_components          = n_components,
    weights               = rep(1 / n_components, n_components),
    mm                    = build_emission(measurement_engine, n_components = n_components, ...),
    sm                    = if (!is.null(structural))
      build_emission(structural, n_components = n_components, ...)
    else NULL,
    n_steps               = n_steps,
    correction            = correction,
    # Which rule turned the step-1 posteriors into the assigned-class variable
    # the correction inverts. Kept on the fit so a saved object still says which
    # rule produced it.
    assignment            = assignment,
    sample_weights        = weights,
    weight_type           = weight_type,
    n_eff                 = n_eff,
    strata                = if (is.null(strata)) rep(1L, n_samples) else strata,
    cluster               = if (is.null(cluster)) seq_len(n_samples) else cluster,
    has_survey_design     = has_survey_design,
    # Retain the indicator matrix so plot() can scale continuous indicators
    # against their observed range (copy-on-write keeps this cheap).
    data                  = X,
    # Store the original descriptor so bootstrap.R can re-fit replicates
    # using the same measurement specification. Missing-data resolution is
    # re-applied per replicate, so the stored value is the user's spec, not the
    # resolved "*_nan" form.
    measurement_descriptor = measurement,
    # Record where and how missing data were handled (NA-free fits store a
    # summary with any_missing = FALSE).
    missing_data           = missing_data,
    # Prior strengths (see R/bayes_constants.R). Stored on the fit so that
    # print()/summary() can report a non-default choice and bootstrap replicates
    # inherit it.
    bayes_constants        = bayes_constants
  )
  class(model_state) <- "mixture_model"

  # Push the constants down onto the emissions, recursively, so that every
  # M-step and refine_lbfgs() read one object rather than each holding its own
  # copy of a default.
  model_state$mm <- .attach_bayes_constants(model_state$mm, bayes_constants)
  model_state$sm <- .attach_bayes_constants(model_state$sm, bayes_constants)

  # Mirror the survey design onto the structural sub-model (see
  # .mirror_design_onto_sm in R/stepwise.R).
  model_state <- .mirror_design_onto_sm(model_state)

  if (n_steps == 1) {
    model_state <- fit_em(model_state, X, Y, n_init, max_iter, random_state,
                          refine = refine, warm_start = warm_start)

  } else {
    model_state <- fit_em(model_state, X, NULL, n_init, max_iter, random_state,
                          refine = refine, warm_start = warm_start)

    # Step 1 metrics (measurement model only)
    if (n_steps == 3)
      model_state$step1_metrics <- .step1_metrics(model_state)

    model_state <- .apply_structural_steps(model_state, X, Y, n_steps,
                                           correction, max_iter, se,
                                           assignment = assignment)
  }

  # Class sorting, display names, and combined-model metrics (see
  # .finalize_model_state in R/stepwise.R).
  model_state <- .finalize_model_state(model_state, X, Y, order_by_size)

  # Collapsed-variance check, after sorting so the class numbers it reports are
  # the ones the user will see. See R/gaussian_boundary.R.
  model_state <- .check_gaussian_degeneracy(model_state, X)

  return(model_state)
}

# ==============================================================================
# User-facing front-end for fit_mixture()
# ==============================================================================

# Measurement families whose complete-data descriptor has a missing-data (FIML)
# counterpart, mapping each base descriptor to the variant that masks NA during
# estimation. Descriptors absent from this table (e.g. structural families) have
# no missing-data variant and pass through resolution unchanged.
.nan_variant <- c(
  bernoulli     = "bernoulli_nan",
  binary        = "binary_nan",
  multinoulli   = "multinoulli_nan",
  categorical   = "categorical_nan",
  gaussian_diag = "gaussian_diag_nan",
  continuous    = "continuous_nan",
  gaussian_unit = "gaussian_unit_nan",
  gaussian      = "gaussian_nan",
  poisson       = "poisson_nan",
  count         = "count_nan"
)

# Resolve a measurement descriptor against its data so the estimator matches the
# data: complete columns keep the fast complete-data form, columns containing NA
# switch to the FIML variant that masks missing cells. Descriptors that already
# name a missing-data variant, or that have no variant, are returned unchanged.
#
# For a nested (mixed) measurement model each block is resolved against the
# columns it governs. Blocks are stored in order with consecutive column counts
# (.normalize_measurement() groups them this way), so a running offset maps each
# block to its columns.
.resolve_emission_descriptor <- function(descriptor, X) {
  nan_strings <- unname(.nan_variant)

  upgrade_one <- function(model, cols_have_na) {
    if (!cols_have_na) return(model)
    if (model %in% nan_strings) return(model)          # already a _nan variant
    variant <- .nan_variant[model]
    if (is.na(variant)) return(model)                  # no missing-data variant
    unname(variant)
  }

  if (is.character(descriptor) && length(descriptor) == 1L)
    return(upgrade_one(descriptor, anyNA(X)))

  if (is.list(descriptor)) {
    offset <- 0L
    for (name in names(descriptor)) {
      n_cols <- descriptor[[name]]$n_columns
      cols   <- seq.int(offset + 1L, offset + n_cols)
      descriptor[[name]]$model <-
        upgrade_one(descriptor[[name]]$model,
                    anyNA(X[, cols, drop = FALSE]))
      offset <- offset + n_cols
    }
    return(descriptor)
  }

  descriptor
}

# Translate a user measurement specification into the descriptor the fitting
# engine consumes, returning the (possibly re-grouped) indicator matrix.
#
# Accepts either a single type string (every indicator shares that type) or a
# named list/vector mapping a measurement type to the indicators it governs,
# for mixed-type models. Indicators may be referenced by column name or index:
#
#   measurement = "binary"
#   measurement = list(binary = c("q1", "q2"), continuous = c("score1"))
#   measurement = list(binary = 1:3, continuous = 4:5)
#
# The suggestion in the missing-argument error, and nowhere else. This is a hint
# for the message text; it must never be used to choose a model. That is the
# whole point of requiring the argument: the storage mode of a column does not
# determine its measurement model.
.measurement_hint <- function(X) {
  if (is.null(X)) return(NULL)
  M <- try(as.matrix(X), silent = TRUE)
  if (inherits(M, "try-error")) return(NULL)
  v <- suppressWarnings(as.numeric(M))
  v <- v[is.finite(v)]
  if (!length(v)) return(NULL)
  ncat <- apply(M, 2, function(z) length(unique(z[!is.na(z)])))
  if (all(ncat == 2L))                                            "binary"
  else if (all(v == round(v)) && all(v >= 0) && all(ncat <= 10L))  "categorical"
  else if (all(v == round(v)) && all(v >= 0))                      "count"
  else                                                             "continuous"
}

# Raised when `measurement` is not supplied. There is deliberately no default:
# a 1-5 column is a legitimate categorical, continuous or count indicator and
# the class solution differs across the three, so a default would settle a
# modelling question by inspecting storage mode, and would make a script's
# meaning depend on the data it happens to be run against.
.require_measurement <- function(X) {
  hint <- .measurement_hint(X)
  msg  <- paste0(
    "`measurement` must be specified. Valid types: \"binary\", ",
    "\"categorical\", \"continuous\", \"count\".")
  if (!is.null(hint)) {
    M    <- as.matrix(X)
    ncat <- apply(M, 2, function(z) length(unique(z[!is.na(z)])))
    why  <- switch(
      hint,
      binary      = "all take two values",
      categorical = sprintf("all take at most %d distinct whole-number values",
                            max(ncat)),
      count       = "all are non-negative whole numbers",
      continuous  = "are not all whole numbers")
    msg <- paste0(msg, sprintf(
      "\nYour %d indicator columns %s, so you probably want\n  measurement = \"%s\"",
      ncol(M), why, hint))
  }
  stop(paste0(msg,
              "\nFor mixed types:\n",
              "  measurement = list(binary = 1:5, continuous = 6:8)"),
       call. = FALSE)
}

# Columns are grouped in the order given so the engine's block structure lines
# up; column names are preserved for display.
.normalize_measurement <- function(measurement, indicators) {
  # Level labels are read before data.matrix() destroys them: it maps a factor
  # or character column to 1-based integer codes, so by the next line there is no
  # way to tell the user which of their categories became the "1" whose
  # probability gets reported.
  labels <- .indicator_level_labels(indicators)

  indicators <- if (is.data.frame(indicators)) data.matrix(indicators)
  else as.matrix(indicators)

  if (is.character(measurement) && length(measurement) == 1L) {
    rec <- .recode_binary(measurement, indicators, seq_len(ncol(indicators)),
                          labels)
    return(list(descriptor = measurement, indicators = rec$indicators,
                recode = rec$map))
  }

  if (!is.list(measurement))
    stop("`measurement` must be a single type string or a named list mapping ",
         "measurement types to indicator columns.", call. = FALSE)

  parts <- as.list(measurement)
  if (is.null(names(parts)) || any(names(parts) == ""))
    stop("For a mixed measurement model, `measurement` must be a named list ",
         "whose names are measurement types (e.g. \"binary\", \"continuous\").",
         call. = FALSE)

  col_names <- colnames(indicators)
  resolve_cols <- function(sel) {
    if (is.character(sel)) {
      if (is.null(col_names))
        stop("Indicator columns were referenced by name, but `indicators` has ",
             "no column names.", call. = FALSE)
      idx <- match(sel, col_names)
      if (anyNA(idx))
        stop("Unknown indicator column name(s): ",
             paste(sel[is.na(idx)], collapse = ", "), call. = FALSE)
      idx
    } else {
      idx <- as.integer(sel)
      if (anyNA(idx) || any(idx < 1L) || any(idx > ncol(indicators)))
        stop("Indicator column indices in `measurement` are out of range.",
             call. = FALSE)
      idx
    }
  }

  keys        <- make.unique(names(parts), sep = "_")
  ordered_idx <- integer(0)
  descriptor  <- list()
  recode      <- list()
  for (i in seq_along(parts)) {
    cols <- resolve_cols(parts[[i]])
    if (any(cols %in% ordered_idx))
      stop("An indicator column was assigned to more than one measurement type.",
           call. = FALSE)
    # Recoding is applied per block, so a mixed model's binary items are covered
    # too. The single-string branch above used to be the only path with any 0/1
    # handling at all, which left a 1/2-coded item in a list spec reaching the
    # Bernoulli likelihood and producing a silently wrong answer.
    rec <- .recode_binary(names(parts)[i], indicators, cols, labels)
    indicators <- rec$indicators
    recode     <- c(recode, rec$map)
    ordered_idx        <- c(ordered_idx, cols)
    descriptor[[keys[i]]] <- list(model = names(parts)[i], n_columns = length(cols))
  }

  unassigned <- setdiff(seq_len(ncol(indicators)), ordered_idx)
  if (length(unassigned) > 0)
    stop("Every indicator column must be assigned a measurement type. ",
         "Unassigned column(s): ", paste(unassigned, collapse = ", "),
         call. = FALSE)

  list(descriptor = descriptor,
       indicators  = indicators[, ordered_idx, drop = FALSE],
       recode      = recode)
}

# Level labels for every column of the user's indicators, so a recoded item can
# say which of the original categories is the one whose probability is reported.
# NULL for a numeric column, whose values speak for themselves.
.indicator_level_labels <- function(indicators) {
  if (!is.data.frame(indicators)) {
    if (is.character(indicators) || is.factor(indicators))
      indicators <- as.data.frame(indicators, stringsAsFactors = FALSE)
    else return(NULL)
  }
  lab <- lapply(indicators, function(v) {
    if (is.factor(v))    return(levels(v))
    if (is.character(v)) return(sort(unique(v[!is.na(v)])))
    if (is.logical(v))   return(c("FALSE", "TRUE"))
    NULL
  })
  names(lab) <- names(indicators)
  lab
}

.binary_families <- c("binary", "bernoulli", "binary_nan", "bernoulli_nan")

# Put binary indicators on the 0/1 scale the Bernoulli emission requires.
#
# The emission is arithmetic, not a lookup -- it computes x*log(p) +
# (1-x)*log(1-p) -- so it needs literal zeros and ones. Everything else an
# applied researcher plausibly has is accepted and converted here: a two-level
# factor or character (first level 0, second 1), a logical, and any numeric pair
# whatever its values, the lower becoming 0. Data already in {0, 1} is untouched
# and silent. Three or more distinct values are left alone for the validator
# downstream to reject, since the right advice there is to use "categorical",
# not to guess at a dichotomy.
#
# The mapping is returned as well as applied. Which level became the 1 decides
# what every reported probability means, and after data.matrix() there is no way
# for the user to recover it from the fit.
.recode_binary <- function(family, indicators, cols, labels = NULL) {
  if (!is.character(family) || length(family) != 1L ||
      !family %in% .binary_families)
    return(list(indicators = indicators, map = list()))

  nms <- colnames(indicators)
  map <- list()
  for (j in cols) {
    v   <- indicators[, j]
    obs <- v[!is.na(v)]
    if (!length(obs)) next
    vals <- sort(unique(obs))
    if (length(vals) != 2L || all(vals == c(0, 1))) next

    nm  <- nms[j] %||% paste0("column ", j)
    lev <- if (!is.null(labels)) labels[[nm]] else NULL
    # data.matrix() codes a factor's levels 1..M in order, so the two observed
    # codes index straight into the recorded labels.
    shown <- if (!is.null(lev) && all(vals %in% seq_along(lev))) lev[vals]
             else format(vals)

    indicators[, j] <- as.numeric(v == vals[2L])
    map[[nm]] <- list(item = nm, zero = shown[1L], one = shown[2L])
  }

  if (length(map)) {
    pairs <- vapply(map, function(m) sprintf("%s (%s -> 0, %s -> 1)",
                                             m$item, m$zero, m$one), character(1))
    message("Recoded binary indicators to 0/1: ",
            paste(pairs, collapse = "; "),
            ". Reported probabilities are of the value shown as 1.")
  }
  list(indicators = indicators, map = map)
}

# Decide whether a distal outcome is continuous or categorical when the user
# leaves outcome_type = "auto". Factors, characters, and integer-valued numerics
# with few distinct values are treated as categorical.
.resolve_outcome_type <- function(outcome, outcome_type) {
  if (outcome_type != "auto") return(outcome_type)
  if (is.factor(outcome) || is.character(outcome)) return("categorical")
  v <- as.numeric(outcome)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return("continuous")
  if (length(unique(v)) <= 10L && all(abs(v - round(v)) < 1e-8))
    return("categorical")
  "continuous"
}

# Pick a display label for the outcome column.
.outcome_label <- function(outcome) {
  if (!is.null(dim(outcome)) && !is.null(colnames(outcome)))
    return(colnames(outcome)[1])
  "outcome"
}

# Best-guess a variable name from the expression a user wrote, covering the
# common ways a single column is referenced: a bare symbol (age), extraction
# with `$` (data$age) or `[[` (data[["age"]]), and single-bracket indexing with
# a character column (data[, "age"]). Returns NULL when no name can be read off
# the expression (e.g. a positional index or a computed vector).
.derive_name <- function(expr) {
  if (is.symbol(expr)) return(as.character(expr))
  if (is.call(expr)) {
    op <- as.character(expr[[1]])
    if (op == "$")  return(as.character(expr[[3]]))
    if (op == "[[" && is.character(expr[[3]])) return(expr[[3]])
    if (op == "[") {
      args <- as.list(expr)
      for (i in seq_along(args)) {
        if (i <= 2L) next                 # skip `[` and the indexed object
        if (is.character(args[[i]])) return(args[[i]])
      }
    }
  }
  NULL
}

# Normalise a predictors / outcome_covariates argument so that prepare_covariates
# always receives something with proper column names. Matrices and data frames
# pass through unchanged (their names are kept); a single bare vector or factor
# is wrapped into a one-column data frame named from the user's expression, with
# `fallback` used only when no name can be inferred.
.as_named_covariates <- function(value, expr, fallback) {
  if (is.null(value) || !is.null(dim(value))) return(value)
  nm  <- .derive_name(expr)
  if (is.null(nm) || !nzchar(nm)) nm <- fallback
  out <- data.frame(value, check.names = FALSE, stringsAsFactors = FALSE)
  names(out) <- nm
  out
}

#' Fit a Latent Class or Latent Profile Mixture Model
#'
#' @description
#' Fits a finite mixture (latent class / latent profile) model. The latent
#' classes are defined by a set of measurement \code{indicators}; optionally, a
#' structural model relates the classes to external variables — either as
#' \code{predictors} of class membership, or as a distal \code{outcome} caused
#' by the classes.
#'
#' @details
#' \strong{The EM convergence rule is fixed and not user-adjustable.} Each
#' random start stops once the log-likelihood improves by less than
#' \code{1e-4} in absolute terms or a relative \code{1e-8}, whichever is
#' looser, and there is no argument that changes this. A looser rule was
#' tried and measured to cost real log-likelihood — on a validated 4-class
#' binary example, stopping early and then polishing with a numerical
#' optimizer recovered only a ninth of what a tighter EM rule found outright
#' — and, more importantly, it degrades the multi-start search itself: a
#' loose rule ranks candidate starts on log-likelihood values too coarse to
#' tell a good basin from a mediocre one, so \code{n_init} stops doing its
#' job. If a fit seems slow, raise \code{n_init} or trim \code{max_iter}
#' rather than looking for a tolerance argument; there is not one to find.
#'
#' @param indicators Matrix or data frame of measurement items that define the
#'   latent classes (rows are observations, columns are items).
#' @param n_classes Number of latent classes/profiles to estimate.
#' @param measurement What kind of variables the indicators are. Either one type
#'   string for all of them — \code{"binary"}, \code{"categorical"},
#'   \code{"continuous"}, \code{"gaussian"}, \code{"count"} — or, for a
#'   mixed-type model, a named list mapping each type to the columns it governs
#'   by name or index, e.g.
#'   \code{list(binary = c("q1","q2"), continuous = "score")}.
#'
#'   \strong{Required; there is no default.} The storage mode of a column does
#'   not determine its measurement model: a 1-5 column is a legitimate
#'   \code{"categorical"}, \code{"continuous"} or \code{"count"} indicator, and
#'   the class solution differs across the three. Inferring the type from the
#'   data would settle a modelling question by inspecting storage mode, and
#'   would make a script's meaning depend on the data it happens to be run
#'   against; a constant default is the same guess with the data-dependence
#'   removed. Omitting the argument is an error that lists the valid types and
#'   suggests one based on your columns — as a hint for you to confirm, not a
#'   choice the package makes.
#'
#'   \strong{Missing values need no special handling.} Any indicator containing
#'   \code{NA} is estimated by full-information maximum likelihood under the
#'   usual missing-at-random assumption, chosen automatically; cases missing on
#'   every indicator are dropped, with a warning saying how many. There is
#'   nothing to switch on. (The \code{"*_nan"} spellings are still accepted, and
#'   force the same estimator, but they are a leftover and you do not need them.)
#'
#'   \strong{Binary items need not be coded 0/1.} A two-level factor or
#'   character, a logical, or any pair of numbers — 1/2 is the commonest — is
#'   converted for you, with a message saying which level became the 1, since
#'   that is what the reported probabilities refer to.
#'
#'   \code{"categorical"} expects integer codes running 1, 2, 3, ...; add 1 to a
#'   0-based item. \code{"count"} fits a Poisson model, one rate per item and
#'   class, and needs non-negative integers.
#' @param predictors Optional covariates that predict latent class membership.
#'   Supplying this fits a class-membership regression (the "predict class"
#'   structural model). Mutually exclusive with \code{outcome}.
#' @param outcome Optional distal outcome caused by the latent classes.
#'   Mutually exclusive with \code{predictors}.
#' @param outcome_covariates Optional covariates that adjust the distal
#'   \code{outcome}.
#' @param outcome_type One of \code{"auto"}, \code{"continuous"}, or
#'   \code{"categorical"}. With \code{"auto"} (default) the type is inferred
#'   from \code{outcome}.
#' @param slopes When \code{outcome_covariates} are supplied, whether their
#'   effect is \code{"pooled"} (one slope shared across classes) or
#'   \code{"class_specific"} (a separate slope per class).
#' @param group Optional observed grouping variable for a multiple-group
#'   model (Collins & Lanza, 2010, sec. 5.7-5.12), e.g. grade or gender.
#'   Unlike \code{predictors}, which only ever shifts class membership, a
#'   group can also be allowed to shift the item-response probabilities
#'   themselves; see \code{group_effects}.
#'
#'   \strong{A covariate whose effect on class membership differs by group}
#'   -- moderation, not just adjustment -- does not need \code{group} at all:
#'   build the interaction into \code{predictors} directly, e.g.
#'   \code{predictors = model.matrix(~ grade * factor(year))[, -1]} handed to
#'   \code{fit_mixture()} in place of a plain covariate matrix. Leave
#'   \code{group} unset for this (or use
#'   \code{group_effects = "measurement"} if the item parameters should also
#'   be free by group): a \code{"prevalence"} or \code{"both"} group effect
#'   appends the group's own design to \code{predictors} internally, so the
#'   group columns would then appear twice. A covariate matrix built this way
#'   also carries none of the column-to-variable bookkeeping a data-frame
#'   \code{predictors} does, so functions like \code{\link{wald_test}} need
#'   the individual interaction column names rather than a single variable
#'   name.
#' @param group_effects What you are allowing the groups to differ in. Pick by
#'   the question you are asking:
#'
#'   \describe{
#'     \item{\code{"prevalence"}}{Do the groups differ in \emph{how many} people
#'       fall in each class? The classes themselves mean the same thing in every
#'       group; only their sizes move. This is the model most applied analyses
#'       want.}
#'     \item{\code{"both"} (default)}{Does each group need its \emph{own} class
#'       solution? Both the class sizes and what the classes look like are free
#'       to differ — the configural model.}
#'     \item{\code{"none"}}{Ignore the grouping variable when estimating.}
#'   }
#'
#'   \strong{The usual workflow} is to fit \code{"prevalence"} and \code{"both"}
#'   and compare them with [`lr_test()`]. That is the test of measurement
#'   invariance (Collins & Lanza, 2010, sec. 5.8): it asks whether letting the
#'   classes mean different things in different groups buys a significantly
#'   better fit. If it does not — the common outcome, and the one you want — keep
#'   \code{"prevalence"} and report the group differences in class sizes, which
#'   are then comparable across groups because the classes are. Comparing
#'   \code{"none"} against \code{"prevalence"} tests whether the sizes differ at
#'   all (sec. 5.11). Pass \code{n_steps = 1} for both fits.
#'
#'   There is a fourth setting, \code{"measurement"}, which frees the item
#'   parameters while holding the class sizes pooled. It answers an unusual
#'   question and is rarely what is wanted; \code{"both"} is the configural model
#'   the invariance literature actually compares against.
#' @section Multiple-group models:
#' Under \code{group_effects = "both"} or \code{"measurement"} each group's item
#' parameters are estimated from that group's own cases, so the \emph{class
#' labels} need not line up across groups: "Class 1" in one group's profile is
#' not guaranteed to be the same kind of class as "Class 1" in another's. Read
#' the per-group profiles by their item patterns rather than by position. The
#' likelihood-ratio comparison is unaffected, since it only differences total
#' log-likelihoods.
#'
#' The log-likelihood is that of the indicators \emph{given} the group; the
#' grouping variable's own distribution is not modelled and its proportions are
#' not counted as parameters. Software that treats a known grouping variable as
#' a latent class variable observed without error adds both. For
#' comparison with such output, a \code{group} fit also carries
#' \code{metrics$ll_knownclass} and \code{metrics$n_params_knownclass}; the
#' difference is a fixed constant and cancels in [`lr_test()`].
#' @param group_invariant_items Item indices or names held equal across
#'   groups even when \code{group_effects} frees the measurement model
#'   (Collins & Lanza's partial-invariance models, sec. 5.9). \code{NULL}
#'   (the default) leaves every item free, i.e. a fully configural model.
#' @param group_invariant_params For continuous indicators, which \emph{kind}
#'   of parameter is held equal across groups when \code{group_effects} frees
#'   the measurement model: \code{"means"}, \code{"covariances"}, or both.
#'   \code{NULL} (the default) frees everything. Where
#'   \code{group_invariant_items} shares whole items, this shares one parameter
#'   matrix for every item — which is the model the latent-profile invariance
#'   literature actually fits: Olivera-Aguilar and Rikoon's (2018)
#'   "unconstrained" model frees the class means across groups while holding the
#'   indicator variances invariant, i.e.
#'   \code{group_invariant_params = "covariances"}, and it is that model, not
#'   the fully heterogeneous one, that their invariance test compares against.
#'   The two are alternatives, not combinable. Categorical indicators have a
#'   single kind of parameter, so for them this constraint and
#'   \code{group_invariant_items} coincide and only the latter is offered.
#' @param variances_equal Logical, for continuous indicators only: hold each
#'   item's variance equal across the classes, so the classes differ in location
#'   only. This is the homoscedastic latent profile model and the default
#'   parameterisation of several commercial programs, and it combines with
#'   \code{group_invariant_params} to give a variance that is free across groups
#'   but shared by the classes within each. Passed through to the measurement
#'   model, so it is also available on an ordinary single-group fit.
#'
#'   \strong{This is the default when \code{measurement} is continuous}, and
#'   \code{variances_equal = FALSE} recovers the class-varying parameterisation
#'   the package used previously. The reason is that the unrestricted
#'   normal-mixture likelihood is \emph{unbounded}: send a class mean to any
#'   single data point and that class's variance to zero and the likelihood
#'   diverges, so no maximum likelihood estimate exists and what the EM reports
#'   is a local optimum whose properties are not guaranteed. Holding the
#'   variances equal across classes bounds the likelihood, and the constrained
#'   estimator is consistent (Day, 1969; Hathaway, 1985). Freeing them also
#'   invites classes that describe non-normality in a single population rather
#'   than distinct subgroups (Bauer & Curran, 2003).
#'
#'   The restriction is substantive and testable, and you are expected to fit
#'   both and compare rather than accept either blindly. It is not offered as
#'   the \emph{safe} choice but as the well-posed one: when the homoscedastic
#'   model is wrong, it fails visibly — a genuinely heteroscedastic class splits
#'   into two, and the comparison against the free model says so. When the free
#'   model is wrong it fails silently, as a boundary solution that gets written
#'   up as a finding.
#' @param n_steps Estimation strategy: 1 (simultaneous), 2, or 3 (recommended
#'   when a structural model is present). Defaults to 3 when \code{predictors}
#'   or \code{outcome} is supplied and left unset, otherwise 1.
#' @param correction Bias correction for 3-step estimation: \code{"none"},
#'   \code{"ML"}, or \code{"BCH"}. When left unset for a 3-step structural
#'   model, a recommended default is chosen (ML for predictors and categorical
#'   outcomes, BCH for continuous outcomes).
#' @param assignment How step 1's posteriors are turned into the assigned-class
#'   variable whose classification error the correction inverts:
#'   \code{"proportional"} (default; each case carries its posterior probability
#'   in every class) or \code{"modal"} (each case is assigned to its most likely
#'   class outright). See \code{\link{add_covariates}}.
#' @param weights,strata,cluster Optional survey design: sampling
#'   \code{weights}, and \code{strata}/\code{cluster} identifiers enabling
#'   design-based (linearization) standard errors.
#' @param weight_type What the numbers in \code{weights} mean.
#'   \code{"sampling"} (the default) treats them as survey or probability
#'   weights, saying how much of the population each case represents; only their
#'   relative sizes matter, so they are rescaled to sum to the number of cases.
#'   \code{"frequency"} treats them as counts of identical cases, as in a
#'   response-pattern table where one row stands for many respondents; the
#'   sample size behind AIC and BIC is then the sum of the counts. Getting this
#'   wrong changes BIC but not the parameter estimates.
#' @param n_init,max_iter,random_state,order_by_size,refine Estimation
#'   controls: number of random starts (default 20), maximum EM iterations, RNG
#'   seed, whether to order classes by size, and whether to run L-BFGS
#'   refinement.
#'
#'   A mixture likelihood usually has several local maxima, so a single start is
#'   a coin toss rather than an estimate; the fit reports how many of the starts
#'   reached the solution it kept. If that count is 1, refit with
#'   \code{n_init = 100} — a maximum found once may simply be the best of a
#'   small sample of the surface, and if it does not replicate at 100 starts the
#'   problem is more likely the specification than the search (Hipp & Bauer,
#'   2006). Set \code{random_state} to make the search reproducible.
#'
#'   The default of 20 is a floor for small, well-separated models rather than a
#'   guarantee: under a correct specification the maximum is typically found by
#'   39\% to 77\% of the starts, so 20 is enough to see a stable maximum when
#'   there is one, but local optima multiply with more classes and poorer
#'   separation. \code{max_iter} defaults to 1000, following the iteration
#'   budget in Biernacki et al. (2003); when a fit stops at the cap, double it.
#'   \code{vignette("estimation")} gives the figures behind both numbers.
#'
#'   The cost is roughly linear in \code{n_init}, so it is easy to budget for a
#'   larger search: \code{n_init = 20} took about 16 seconds on a 4-class,
#'   7-item binary model with 2,587 cases, so \code{n_init = 200} on the same
#'   model is a couple of minutes and \code{n_init = 1000} is closer to ten.
#'
#'   More starts are not \emph{monotonically} better on a large model. Where each
#'   EM iteration is expensive the search runs in two stages: a short pass ranks
#'   the starts and only the best fraction is run to convergence. That ranking is
#'   taken before the starts have converged, so it can discard a slow-climbing
#'   basin that would have won, and raising \code{n_init} changes which starts
#'   survive rather than simply adding to them. A larger \code{n_init} can
#'   therefore land on a slightly worse solution than a smaller one. On a big
#'   model, compare two or three values rather than assuming the largest is best.
#' @param bayes_constants Optional list adjusting the strength of the weak
#'   priors the estimator places on each block of parameters. Named
#'   \code{latent} (class weights and, in the transition models, the initial
#'   and transition probabilities), \code{categorical} (item-response
#'   probabilities), \code{poisson} (count rates), and \code{variances}
#'   (class-specific variances of continuous indicators); all default to
#'   \code{1}. Each is a number of pseudo-observations spread over the classes,
#'   so its influence shrinks as the sample grows.
#'
#'   The default of one observation, and the fact that it is spread in agreement
#'   with each item's own observed marginal rather than uniformly over the
#'   cells, are both taken from Galindo Garre and Vermunt (2006). Their
#'   simulation compares Jeffreys, normal and two Dirichlet priors against
#'   maximum likelihood and the parametric bootstrap, at a prior strength of a
#'   single added case, and the marginal-preserving Dirichlet has the lowest
#'   root median squared error and the best coverage in every condition they
#'   study — most clearly where a true parameter is near the boundary and small
#'   samples send the maximum likelihood estimate to infinity. The constant
#'   Dirichlet, which spreads the same one observation uniformly, degrades as
#'   the table grows: with nine binary items that observation is spread over 512
#'   cells rather than 32, and its estimates shrink too far toward zero. This is
#'   why \code{categorical} is not simply an add-one adjustment.
#'
#'   The \code{variances} prior is centred on each item's own observed marginal
#'   variance rather than on a fixed number, which is what makes it mean the same
#'   thing on a five-point scale and on an income variable: rescaling an
#'   indicator rescales the prior with it. That invariance is the reason the
#'   penalised-likelihood literature recommends a data-scaled penalty over a
#'   constant one (Chen, Tan, & Zhang, 2008, sec. 4).
#'
#'   Two uses. \strong{Reproducing an unregularized fit:} setting a constant to
#'   \code{0} removes that prior and gives plain maximum likelihood for that
#'   block. This is an escape hatch for matching a reference analysis, not a
#'   recommended setting — the unpenalised mixture likelihood for a mixture of
#'   normals is unbounded, so a \emph{global} maximum likelihood estimate does
#'   not exist (Day, 1969; Kiefer & Wolfowitz, 1956) and what an unregularized
#'   program reports is a local maximum.
#'
#'   \strong{Rescuing a collapsed fit:} this situation is now rare, because a
#'   continuous measurement model holds the variances equal across classes by
#'   default and that bounds the likelihood; it arises when you have explicitly
#'   passed \code{variances_equal = FALSE}. If a fit warns that a class variance
#'   has collapsed, the prior is one of three remedies, and not the first to
#'   reach for. Constraining the variances to be equal across classes
#'   (\code{variances_equal = TRUE}) bounds the likelihood so the problem cannot
#'   arise at all; fitting fewer classes often removes the class that was
#'   describing a spike. Where neither is acceptable substantively, raise
#'   \code{variances}. A useful starting point is \strong{one artificial
#'   observation per class}, i.e. \code{variances = n_classes}, increasing it a
#'   little at a time if the warning persists. State it that way rather than as a
#'   bare number: the constant is divided among the classes, so
#'   \code{variances = n_classes} is what holds the prior at one
#'   pseudo-observation per class at \emph{every} number of classes — the way the
#'   penalised-mixture literature applies it, the same amount to every component
#'   — while a fixed value drifts as classes are added. Then check the result
#'   rather than just that the warning stopped — the flagged variance should no
#'   longer be far below the others in the model, and its class mean should have
#'   come off the floor or ceiling of the response scale. It is worth looking at
#'   the distribution of the named item as well, for the floor, ceiling or spike
#'   the class latched onto.
#'
#'   Raising \code{n_init} is \emph{not} a remedy here. A collapsed variance is
#'   not a convergence failure but a property of the likelihood, which really is
#'   unbounded in that direction, so a longer search can find a taller spike and
#'   report a better log-likelihood for a worse solution. For the same reason a
#'   flagged fit's BIC is not comparable with a clean fit's.
#'
#'   The theory behind the prior settles its form and not its size: within the
#'   conditions Chen, Tan, & Zhang (2008) impose — which any fixed positive
#'   constant meets — a penalty that diverges as a variance approaches zero and
#'   grows more slowly than the sample yields a consistent estimator whatever its
#'   constant. That result is proved for \emph{univariate} normal mixtures; the
#'   multivariate case, which is what this package fits, they leave open (sec. 5).
#'   The choice of constant is therefore all the more a finite-sample judgement,
#'   which is why the default is deliberately weak and the guidance above is a
#'   rule with a check rather than a magic value.
#'
#'   This is not a tuning menu. The defaults are the intended settings.
#' @param se How standard errors for \code{predictors} are computed in a 2- or
#'   3-step model. \code{"corrected"} (the default) adds the variance carried
#'   over from step 1 to the step-3 sandwich, following Bakk et al.
#'   (2014); \code{"robust"} reports the sandwich alone;
#'   \code{"hessian"} inverts the step-3 observed
#'   information only. See \code{\link{covariate_se}}.
#' @param X,Y,n_components,structural Deprecated legacy arguments retained for
#'   backward compatibility; prefer \code{indicators}, \code{n_classes},
#'   \code{predictors}, and \code{outcome}.
#' @param ... Passed through to the measurement-model constructors.
#'
#' @return A fitted \code{mixture_model} object.
#'
#' @references
#' Biernacki, C., Celeux, G., & Govaert, G. (2003). Choosing starting values for
#' the EM algorithm for getting the highest likelihood in multivariate Gaussian
#' mixture models. \emph{Computational Statistics & Data Analysis},
#' \emph{41}(3-4), 561-575. \doi{10.1016/S0167-9473(02)00163-9}
#' (the iteration budget behind \code{max_iter}).
#'
#' Hipp, J. R., & Bauer, D. J. (2006). Local solutions in the estimation of
#' growth mixture models. \emph{Psychological Methods}, \emph{11}(1), 36-53.
#' \doi{10.1037/1082-989X.11.1.36} (the replication rates behind
#' \code{n_init}).
#'
#' Shireman, E., Steinley, D., & Brusco, M. J. (2017). Examining the effect of
#' initialization strategies on the performance of Gaussian mixture modeling.
#' \emph{Behavior Research Methods}, \emph{49}(1), 282-293.
#' \doi{10.3758/s13428-015-0697-6}
#'
#' Collins, L. M., & Lanza, S. T. (2010). \emph{Latent Class and Latent
#' Transition Analysis: With Applications in the Social, Behavioral, and Health
#' Sciences}. Wiley.
#'
#' Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. In
#' T. D. Little (Ed.), \emph{The Oxford Handbook of Quantitative Methods}
#' (Vol. 2, pp. 551-611). Oxford University Press.
#'
#' Vermunt, J. K., & Magidson, J. (2002). Latent class cluster analysis. In
#' J. A. Hagenaars & A. L. McCutcheon (Eds.), \emph{Applied Latent Class
#' Analysis} (pp. 89-106). Cambridge University Press.
#'
#' Galindo Garre, F., & Vermunt, J. K. (2006). Avoiding boundary estimates in
#' latent class analysis by Bayesian posterior mode estimation.
#' \emph{Behaviormetrika}, \emph{33}(1), 43-59. \doi{10.2333/bhmk.33.43}
#' (the priors behind \code{bayes_constants}).
#'
#' Chen, J., Tan, X., & Zhang, R. (2008). Inference for normal mixtures in mean
#' and variance. \emph{Statistica Sinica}, \emph{18}(2), 443-465
#' (the penalty behind \code{bayes_constants$variances}).
#'
#' Day, N. E. (1969). Estimating the components of a mixture of normal
#' distributions. \emph{Biometrika}, \emph{56}(3), 463-474
#' (the unbounded likelihood behind the \code{variances_equal} default).
#'
#' Hathaway, R. J. (1985). A constrained formulation of maximum-likelihood
#' estimation for normal mixture distributions. \emph{The Annals of
#' Statistics}, \emph{13}(2), 795-800 (consistency of the constrained
#' estimator).
#'
#' Bauer, D. J., & Curran, P. J. (2003). Distributional assumptions of growth
#' mixture models: Implications for overextraction of latent trajectory
#' classes. \emph{Psychological Methods}, \emph{8}(3), 338-363.
#' \doi{10.1037/1082-989X.8.3.338}
#'
#' Kiefer, J., & Wolfowitz, J. (1956). Consistency of the maximum likelihood
#' estimator in the presence of infinitely many incidental parameters.
#' \emph{The Annals of Mathematical Statistics}, \emph{27}(4), 887-906.
#'
#' McLachlan, G. J., & Peel, D. (2000). \emph{Finite Mixture Models}. Wiley
#' (on the unbounded likelihood and spurious solutions).
#'
#' Olivera-Aguilar, M., & Rikoon, S. H. (2018). Assessing measurement invariance
#' in multiple-group latent profile analysis. \emph{Structural Equation
#' Modeling}, \emph{25}(3), 439-452. \doi{10.1080/10705511.2017.1408015}
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#'
#' \dontrun{
#' # Class membership predicted by a covariate (3-step, ML by default)
#' fit_mixture(X, n_classes = 2, measurement = "binary", predictors = age)
#'
#' # Distal outcome with a class-specific covariate slope
#' fit_mixture(X, n_classes = 2, measurement = "binary", outcome = bmi,
#'             outcome_covariates = age, slopes = "class_specific")
#'
#' # Mixed-type indicators
#' fit_mixture(items, n_classes = 3,
#'             measurement = list(binary = 1:5, continuous = 6:8))
#' }
#'
#' @export
fit_mixture <- function(indicators = NULL,
                        n_classes = 2,
                        measurement,
                        predictors = NULL,
                        outcome = NULL,
                        outcome_covariates = NULL,
                        outcome_type = c("auto", "continuous", "categorical"),
                        slopes = c("pooled", "class_specific"),
                        group = NULL,
                        group_effects = c("both", "measurement", "prevalence", "none"),
                        group_invariant_items = NULL,
                        group_invariant_params = NULL,
                        variances_equal = NULL,
                        n_steps = 1,
                        correction = "none",
                        assignment = c("proportional", "modal"),
                        n_init = 20,
                        max_iter = 1000,
                        random_state = NULL,
                        order_by_size = TRUE,
                        weights = NULL,
                        weight_type = c("sampling", "frequency"),
                        strata = NULL,
                        cluster = NULL,
                        refine = TRUE,
                        bayes_constants = NULL,
                        se = c("corrected", "robust", "hessian"),
                        X = NULL, Y = NULL, n_components = NULL, structural = NULL,
                        ...) {

  se            <- match.arg(se)
  assignment    <- match.arg(assignment)
  outcome_type  <- match.arg(outcome_type)
  slopes        <- match.arg(slopes)
  group_effects <- match.arg(group_effects)
  weight_type   <- match.arg(weight_type)
  steps_set     <- !missing(n_steps)
  corr_set      <- !missing(correction)
  measurement_missing <- missing(measurement)

  # Capture the unevaluated expressions the user supplied so that a single
  # covariate passed as a bare vector (e.g. data$age or data[, "age"]) can be
  # given an informative name; see .as_named_covariates() below.
  predictors_expr <- substitute(predictors)
  outcome_cov_expr <- substitute(outcome_covariates)

  # --- Legacy interface bridge ------------------------------------------------
  legacy <- !is.null(X) || !is.null(Y) || !is.null(n_components) ||
    !is.null(structural)
  if (!is.null(X) && is.null(indicators)) indicators <- X
  if (!is.null(n_components))              n_classes  <- n_components

  # Checked here rather than in the signature: there is no default, because the
  # storage mode of a column does not determine its measurement model. Comes
  # before the legacy bridge below, which would otherwise pass the missing
  # argument straight through to the engine's own "binary" default.
  if (measurement_missing) .require_measurement(indicators)

  if (legacy) {
    message("Note: `X`, `Y`, `n_components`, and `structural` are the legacy ",
            "interface. The current arguments are `indicators`, `n_classes`, ",
            "`predictors`, and `outcome` / `outcome_covariates`.")
    return(fit_mixture_internal(
      X = indicators, Y = Y, n_components = n_classes,
      measurement = measurement, structural = structural,
      n_steps = n_steps, correction = correction, assignment = assignment,
      n_init = n_init,
      max_iter = max_iter, random_state = random_state,
      order_by_size = order_by_size, weights = weights,
      weight_type = weight_type,
      strata = strata, cluster = cluster, refine = refine,
      bayes_constants = bayes_constants, se = se, ...))
  }

  if (is.null(indicators))
    stop("`indicators` is required: the matrix of measurement items that ",
         "define the latent classes.", call. = FALSE)

  if (!is.null(predictors) && !is.null(outcome))
    stop("Specify either `predictors` (to model class membership) or ",
         "`outcome` (a distal outcome), not both in one model. To run both ",
         "analyses from one solution, fit the unconditional model and use ",
         "add_covariates() and add_outcome() on it.", call. = FALSE)
  if (!is.null(outcome_covariates) && is.null(outcome))
    stop("`outcome_covariates` requires an `outcome`.", call. = FALSE)
  if (!is.null(group) && group_effects %in% c("both", "prevalence") &&
      !is.null(outcome))
    stop("`group` with a prevalence effect (`group_effects = \"both\"` or ",
         '"prevalence") uses the same structural-model slot as `outcome`; ',
         "combine `group` with `predictors` instead, or set `group_effects` ",
         'to "measurement" or "none".', call. = FALSE)

  if (!is.null(group_invariant_params) &&
      !group_effects %in% c("both", "measurement"))
    stop("`group_invariant_params` constrains the measurement model across ",
         "groups, which only `group_effects = \"both\"` or \"measurement\" ",
         "frees in the first place.", call. = FALSE)
  # A continuous measurement model defaults to the homoscedastic
  # parameterisation (Day 1969; Hathaway 1985: the unrestricted likelihood is
  # unbounded, so only the constrained problem has an MLE). Resolved here rather
  # than in the signature because the gate below rejects
  # `variances_equal = TRUE` for any other measurement type, and a bare TRUE
  # default would fire it on every categorical fit. An explicit TRUE still
  # reaches the gate and is still rejected; only the unset case is resolved from
  # the descriptor. `gaussian_model()`'s own default deliberately stays FALSE:
  # it is reached by the growth, time-block and group-block paths as well, and
  # flipping it there would silently change LCGA, GMM and RMLCA. Only the
  # user-facing entry points resolve the default.
  if (is.null(variances_equal))
    variances_equal <- is.character(measurement) && length(measurement) == 1L &&
      measurement %in% c("continuous", "continuous_nan",
                         "gaussian_diag", "gaussian_diag_nan")

  if (isTRUE(variances_equal) &&
      !(is.character(measurement) && length(measurement) == 1L &&
        measurement %in% c("continuous", "continuous_nan",
                           "gaussian_diag", "gaussian_diag_nan")))
    stop("`variances_equal` constrains the class variances of continuous ",
         "indicators; it has no meaning for `measurement = ",
         deparse1(measurement), "`. For a mixed measurement model, pass the ",
         "constraint in the measurement descriptor itself.", call. = FALSE)

  # --- Measurement model (single-type or mixed) -------------------------------
  mm                 <- .normalize_measurement(measurement, indicators)
  X_use              <- mm$indicators
  measurement_engine <- mm$descriptor

  # --- Grouping variable: measurement effect (multiple-group model) ----------
  # Reuses the time-blocks trick across a group axis instead of a time axis
  # (R/group_blocks.R): each group gets its own copy of the item block, and a
  # case's other-group blocks are structurally missing, which FIML already
  # handles correctly. The prevalence effect is handled below, alongside
  # `predictors`, since both use the same class-membership regression.
  group_info       <- NULL
  group_extra_args <- list()
  group_warm_start <- NULL
  # The indicator matrix as the user supplied it, kept because
  # .pad_group_blocks() below replaces X_use with a matrix that is two-thirds
  # structurally empty by construction. Everything reported about *missing data*
  # has to describe this one. See .fix_group_block_reporting().
  X_unpadded <- X_use
  if (!is.null(group)) {
    group_info <- .lta_group_design(group, nrow(X_use))

    if (group_effects %in% c("both", "measurement")) {
      # Warm start (see R/group_blocks.R). The pooled measurement model is fitted
      # first, on the unpadded data, and replicated into every group block as one
      # extra restart of the group-varying search. It is the same model the
      # invariance test compares against, so this is what makes the comparison a
      # test rather than a lower bound. Fitted before X_use is padded and
      # `measurement_engine` is overwritten below, since both are needed as they
      # stand here.
      .sub_fit <- function(rows) {
        args <- c(list(X = if (is.null(rows)) X_use else X_use[rows, , drop = FALSE],
                       Y = NULL, n_components = n_classes,
                       measurement = measurement_engine, structural = NULL,
                       n_steps = 1, correction = "none", n_init = n_init,
                       max_iter = max_iter, random_state = random_state,
                       order_by_size = order_by_size,
                       weights = if (is.null(rows)) weights else weights[rows],
                       weight_type = weight_type,
                       strata  = if (is.null(rows)) strata  else strata[rows],
                       cluster = if (is.null(rows)) cluster else cluster[rows],
                       refine = refine, bayes_constants = bayes_constants,
                       se = se, variances_equal = isTRUE(variances_equal)),
                  list(...))
        out <- try(suppressWarnings(suppressMessages(
          do.call(fit_mixture_internal, args))), silent = TRUE)
        if (inherits(out, "try-error")) NULL else out
      }

      pooled_fit <- .sub_fit(NULL)
      # A failure here costs the warm start and nothing else: the group model is
      # still fitted from random starts exactly as before.
      if (!is.null(pooled_fit)) {
        g_idx     <- as.integer(group_info$factor)
        group_mms <- lapply(seq_len(nlevels(group_info$factor)), function(g) {
          rows <- which(g_idx == g)
          # A group too small to identify K classes on its own is left to the
          # pooled parameters rather than fitted badly.
          if (length(rows) < 2L * n_classes) return(NULL)
          one <- .sub_fit(rows)
          if (is.null(one)) NULL else .align_to_pooled(one$mm, pooled_fit$mm)
        })
        group_warm_start <- .group_blocks_warm_start(pooled_fit, group_mms)
      }

      item_names <- colnames(X_use) %||% paste0("Item", seq_len(ncol(X_use)))
      grp_spec <- .resolve_invariance(
        if (is.null(group_invariant_items)) "none" else "partial",
        group_invariant_items, item_names, measurement)
      n_groups <- nlevels(group_info$factor)
      X_grp  <- .pad_group_blocks(X_use, group_info$factor)
      engine <- .longitudinal_measurement_spec(measurement, X_grp,
                                               n_items = ncol(X_use),
                                               n_times = n_groups)
      X_use              <- engine$X
      measurement_engine <- "group_blocks"
      group_extra_args <- list(
        n_items          = length(item_names),
        n_groups         = n_groups,
        sub_model        = engine$sub_model,
        invariant_items  = grp_spec$invariant_items,
        invariant_params = group_invariant_params %||% character(0),
        max_val          = engine$max_val
      )
    }
  }

  # --- Structural model -------------------------------------------------------
  structural_engine <- NULL
  Y_use             <- NULL

  if (!is.null(predictors)) {
    structural_engine <- "predict_class"
    Y_use             <- prepare_covariates(
      .as_named_covariates(predictors, predictors_expr, "predictor"))

  } else if (!is.null(outcome)) {
    outcome_spec      <- .build_outcome_spec(outcome, outcome_covariates,
                                             outcome_type, slopes,
                                             outcome_cov_expr)
    structural_engine <- outcome_spec$engine
    Y_use             <- outcome_spec$Y
  }

  # --- Grouping variable: prevalence effect -----------------------------------
  # Enters `group` as a class-membership covariate, exactly like `predictors`
  # (and combined with it, if both are given): each group gets its own free
  # prevalences while the measurement model above is whatever the
  # `group_effects` "measurement" branch left it (pooled, unless that branch
  # also ran).
  if (!is.null(group) && group_effects %in% c("both", "prevalence")) {
    structural_engine <- "predict_class"
    Y_use <- if (is.null(Y_use)) group_info$design
             else .cbind_covariates(Y_use, group_info$design)
  }

  # --- Friendly defaults when a structural model is present -------------------
  if (!is.null(structural_engine) && !steps_set) {
    n_steps <- 3L
    # `group` reaches this branch through the prevalence effect above, so a user
    # who asked only for a multiple-group measurement model gets a three-step
    # fit without having asked for one. That changes what metrics$ll is, which
    # is easy to miss inside a loop over models, so say so here rather than
    # leaving it to the documentation.
    group_only <- !is.null(group) && is.null(predictors) && is.null(outcome)
    if (group_only) {
      message(paste("Using 3-step estimation (set `n_steps` to override):",
                    "`group` is being modelled as a predictor of class",
                    "membership, so `metrics$ll` is the structural model's",
                    "log-likelihood and is not comparable with a 1-step fit's.",
                    "For a multiple-group measurement model, pass",
                    "`n_steps = 1` explicitly."))
    } else {
      message("Using 3-step estimation (set `n_steps` to override).")
    }
  }
  if (!is.null(structural_engine) && n_steps == 3L && !corr_set) {
    correction <- if (identical(structural_engine, "predict_class")) "ML"
    else if (startsWith(structural_engine, "categorical")) "ML"
    else "BCH"
    message(sprintf("Using '%s' bias correction (set `correction` to override).",
                    correction))
  }

  dots <- c(list(...), list(variances_equal = isTRUE(variances_equal)))
  if (length(group_extra_args)) dots <- utils::modifyList(dots, group_extra_args)

  fit <- do.call(fit_mixture_internal, c(list(
    X = X_use, Y = Y_use, n_components = n_classes,
    measurement = measurement_engine, structural = structural_engine,
    n_steps = n_steps, correction = correction, assignment = assignment,
    n_init = n_init,
    max_iter = max_iter, random_state = random_state,
    order_by_size = order_by_size, weights = weights,
    weight_type = weight_type,
    strata = strata, cluster = cluster, refine = refine,
    bayes_constants = bayes_constants, warm_start = group_warm_start,
    se = se), dots))

  if (length(mm$recode)) fit$binary_recode <- mm$recode

  if (!is.null(group)) {
    fit$group_info    <- group_info
    fit$group_effects <- group_effects
    fit <- .fix_group_block_reporting(fit, X_unpadded, group_info)
    fit <- .add_knownclass_scale(fit, group_info)
  }

  # Non-convergence used to be near-impossible to hit and was reported only by
  # print(), which a user working from summary() or the coefficients never sees.
  # It became reachable when the emissions L-BFGS does not polish were given a
  # stopping rule tight enough to be trusted: those models genuinely need
  # hundreds of iterations, and the default cap is 1000. Silently returning the
  # iterate EM happened to be on at the cap is the one outcome worth
  # interrupting for, so say it here rather than leaving it to be noticed.
  if (isFALSE(fit$converged)) .warn_non_convergence(max_iter)

  # A maximum found by one start out of many is the other outcome worth
  # interrupting for, and it too was reported only by print(). The growth models
  # are the exception: `structured_normal` is fit_gmm()'s call into this
  # wrapper, and at the moment this function returns, fit_gmm() has not yet run
  # its boundary check, so the guard inside .check_replication() would be
  # looking at a flag that is not set yet. fit_gmm() issues the warning itself
  # once it is. A hand-written fit_mixture(measurement = "structured_normal")
  # therefore does not get it; that is the accepted cost of not having two
  # warnings contradict each other.
  if (!identical(measurement, "structured_normal")) .check_replication(fit)

  fit
}

#' Print a Brief Summary of a Fitted Mixture Model
#'
#' @description
#' Prints a compact overview of the fitted model including: number of classes,
#' estimation method, convergence status, the fit indices, and estimated class
#' proportions. The indices shown are the same ones
#' \code{\link{compare_mixtures}} tabulates, so one printed model and a table
#' over a range of K can be read together; lower AIC, BIC and SABIC are better.
#' For full parameter tables, use
#' \code{\link{summary.mixture_model}} (structural parameters) or
#' \code{\link{measurement_summary}} (item parameters).
#'
#' @param x A fitted \code{mixture_model} object returned by
#'   \code{\link{fit_mixture}}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return Invisibly returns \code{x}. Called for its printed side-effect.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' print(fit)
#' # or equivalently:
#' fit
#'
#' @export
print.mixture_model <- function(x, ...) {
  cat("=========================================================\n")
  cat("                  LATENT MIXTURE MODEL                   \n")
  cat("=========================================================\n")
  cat(sprintf("Classes Estimated  : %d\n", x$n_components))
  cat(sprintf("Estimation Method  : %d-step\n", x$n_steps))
  if (x$n_steps == 3)
    cat(sprintf("Correction Method  : %s%s\n", x$correction,
                if (!is.null(x$assignment) && x$correction != "none")
                  sprintf(" (%s assignment)", x$assignment) else ""))
  cat(sprintf("Converged          : %s (in %d iterations)\n", x$converged, x$n_iter))
  # Deleted cases are reported on their own line, ahead of the FIML summary, so
  # the printed sample size can always be reconciled with the input data.
  if (isTRUE(x$missing_data$n_empty_rows > 0L))
    cat(sprintf("Cases Removed      : %d of %d with no observed indicator (n = %d analysed)\n",
                x$missing_data$n_empty_rows, x$missing_data$n_input_rows,
                x$missing_data$n_input_rows - x$missing_data$n_empty_rows))
  # A group-varying measurement model is estimated as one block per group. That
  # is a fact about the model, and is stated as one; it is not missingness, even
  # though the estimator reaches it through the same FIML machinery.
  if (!is.null(x$missing_data$group_blocks)) {
    gb <- x$missing_data$group_blocks
    cat(sprintf("Measurement Model  : group-varying, %d groups x %d items\n",
                gb$n_groups, gb$n_items))
  }
  if (!is.null(x$missing_data) && isTRUE(x$missing_data$any_missing)) {
    md <- x$missing_data
    cat(sprintf("Missing Data       : %d / %d cells (%.1f%%) in %d item%s \u2014 %s\n",
                md$n_missing, md$n_cells, 100 * md$prop_missing,
                md$n_items_affected, if (md$n_items_affected == 1L) "" else "s",
                md$handled_by))
  }
  if (identical(x$weight_type, "frequency"))
    cat(sprintf("Case Weights       : frequency counts (%s cases in %d rows)\n",
                format(x$n_eff), length(x$sample_weights)))
  cat("---------------------------------------------------------\n")
  # A three-step fit carries two sets of metrics. Every line below is read off
  # one of them, never a mixture of the two: pairing a step-1 log-likelihood
  # with a step-3 BIC is worse than printing no criteria at all.
  .print_fit_indices(x$step1_metrics %||% x$metrics,
                     suffix    = if (!is.null(x$step1_metrics)) " (Step 1)" else "",
                     flag_bic  = !is.null(x$degenerate) ||
                       length(x$growth$boundary) > 0L)
  .print_replication_note(x)
  cat("---------------------------------------------------------\n")
  cat("Class Weights (Sizes):\n")
  for (i in seq_along(x$weights))
    cat(sprintf("  Class %d: %.2f%%\n", i, x$weights[i] * 100))
  cat("=========================================================\n")
  # A warning is transient; someone opening a saved fit months later should
  # still see that its variances collapsed. See R/gaussian_boundary.R.
  .print_degenerate_note(x)
  cat("Type summary(model) for structural parameters or measurement_summary(model) for item parameters.\n")
}

#' Compare Mixture Models Across a Range of Class Numbers
#'
#' @description
#' Fits a sequence of measurement-only mixture models, one for each value of
#' \code{k} in \code{k_range}, and returns a table of fit indices to guide
#' class enumeration. The best model according to BIC is identified
#' automatically.
#'
#' @details
#' **Reading the `Entropy` column.** Relative entropy describes how cleanly a
#' solution separates the classes, on a 0-to-1 scale. The usual anchors are 0.40,
#' 0.60 and 0.80 for low, medium and high separation (Clark & Muthen, 2009, as
#' reported by Lee et al., 2023, p. 653); Ram and Grimm (2009, p. 571) put the
#' same point as "high values of entropy (>.80) indicate that individuals are
#' classified with confidence", and suggest preferring the higher-entropy model
#' when choosing among models with similar BIC. Those anchors are on the same
#' normalisation this package uses.
#'
#' Entropy is not evidence for how many classes there are, and there is no
#' threshold it has to clear: "there are no set cut-off criteria for deciding
#' whether the entropy is reasonably high" (Jung & Wickrama, 2008, p. 312). The
#' numbers above are for reading a table, and mixtureEM applies no entropy
#' threshold anywhere.
#'
#' **Reading the `Unreplicated` column.** `TRUE` means the reported maximum for
#' that K was found by exactly one random start. Refit those values of K with
#' more starts (`n_init = 100` is the usual next step) before reporting them;
#' see `vignette("estimation")`.
#'
#' **Reading the `VLMR` columns.** They appear only when `vlmr` is set, and are
#' off by default for two reasons. One is cost: the test needs a numerical
#' Hessian for each model, which is quadratic in the number of parameters, and a
#' function people call casually should not pay that unasked. The other is that
#' the test does not deserve to be printed as a matter of course. Vermunt (2024)
#' concludes that "neither of the two implementations yield uniformly
#' distributed p-values under the correct null hypothesis, indicating this test
#' is not the best model selection tool in mixture modeling".
#'
#' The two implementations differ only in which covariance matrix of the
#' parameters enters the reference distribution: one program uses the ordinary
#' one, another the robust (sandwich) one (Vermunt, 2024). The difference is not
#' cosmetic: on the same data the two can return p = .00 and p = .15. The robust
#' version's reference distribution is much more sensitive to the particular
#' sample, especially when the classes are poorly separated. Neither version's
#' p-values are uniform under the null, so treat a VLMR result as one input
#' among several and prefer [`blrt()`] where it is affordable.
#'
#' The reason is known and is not a numerical artefact. The reference
#' distribution is derived from a theorem (Vuong, 1989, Theorem 3.3) that
#' requires the parameters of the larger model to be identified at the point
#' where it reduces to the smaller one. A mixture never satisfies this: the
#' larger model reproduces the smaller only by emptying a class or by
#' duplicating one, and in each case some parameters vanish from the likelihood
#' and the information matrix is singular (Jeffries, 2003). The test is
#' therefore best read as a descriptive comparison rather than a calibrated
#' p-value.
#'
#' @param X A numeric matrix or data frame of indicator variables.
#' @param k_range Integer vector of class numbers to fit. All values must be >= 1. Default is \code{1:5}.
#' @param measurement Character string or named list specifying the measurement
#'   model type. Required; see \code{\link{fit_mixture}} for the accepted values
#'   and for why there is no default.
#' @param n_init Positive integer. Number of random restarts per model.
#'   Default is \code{10}.
#' @param n_steps Integer. Estimation method: \code{1}, \code{2}, or \code{3}.
#'   Default is \code{1}.
#' @param vlmr Character string. Whether to add the Vuong-Lo-Mendell-Rubin test
#'   of K against K+1 classes, and in which form: \code{"none"} (the default),
#'   \code{"standard"} for Vuong's own formulae on the ordinary covariance
#'   matrix, \code{"robust"} for the variant that substitutes the sandwich
#'   covariance, or \code{"both"}. See Details for why it is off by default.
#' @param ... Additional arguments passed to \code{\link{fit_mixture}}.
#'
#' @return An object of class `mixture_comparison`: a named list with three
#'   elements, which can be indexed exactly as a plain list.
#'   * `fit_table` Data frame with one row per K and columns `Classes`, `LL`,
#'     `Params`, `AIC`, `BIC`, `SABIC`, `Entropy` and `Unreplicated`. With
#'     `vlmr` set it also carries `VLMR_LR` and one p-value column per
#'     requested form (`VLMR_p`, `VLMR_p_robust`); each row tests its own K
#'     against the next one in the table, so the last row is `NA`.
#'   * `models` Named list of fitted `mixture_model` objects, one per K
#'     (names are `"K1"`, `"K2"`, etc.).
#'   * `best_k` Integer. The value of K with the lowest BIC.
#'   * `vlmr` Present only when `vlmr` is set: one entry per row of the table
#'     holding the likelihood-ratio statistic and, for each requested form, the
#'     mean and standard deviation of the reference distribution alongside the
#'     p-value. Those two moments are what say which distribution produced a
#'     given p-value, and are not printed.
#'
#'   [`plot()`][plot.mixture_comparison] draws the criteria against K.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
#'                            n_init = 5)
#' result$fit_table
#' result$best_k
#'
#' @references
#' Nylund, K. L., Asparouhov, T., & Muthen, B. O. (2007). Deciding on the number
#' of classes in latent class analysis and growth mixture modeling: A Monte
#' Carlo simulation study. \emph{Structural Equation Modeling}, \emph{14}(4),
#' 535-569. \doi{10.1080/10705510701575396}
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
#'
#' Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. In
#' T. D. Little (Ed.), \emph{The Oxford Handbook of Quantitative Methods}
#' (Vol. 2, pp. 551-611). Oxford University Press.
#'
#' Vermunt, J. K. (2024). The Vuong-Lo-Mendell-Rubin test for latent class and
#' latent profile analysis. \emph{Methodology}, \emph{20}(1), e12467.
#' \doi{10.5964/meth.12467}
#'
#' Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, \emph{57}(2), 307-333.
#'
#' Lo, Y., Mendell, N. R., & Rubin, D. B. (2001). Testing the number of
#' components in a normal mixture. \emph{Biometrika}, \emph{88}(3), 767-778.
#'
#' Jeffries, N. O. (2003). A note on "Testing the number of components in a
#' normal mixture". \emph{Biometrika}, \emph{90}(4), 991-994.
#'
#' Imhof, J. P. (1961). Computing the distribution of quadratic forms in normal
#' variables. \emph{Biometrika}, \emph{48}(3/4), 419-426.
#'
#' @export
compare_mixtures <- function(X, k_range = 1:5, measurement,
                             n_init = 10, n_steps = 1,
                             vlmr = c("none", "standard", "robust", "both"),
                             ...) {
  if (missing(measurement)) .require_measurement(X)
  vlmr <- match.arg(vlmr)
  # k=0 would silently fit a degenerate model with LL=-Inf; negative
  # k values crash deep in initialisation with a cryptic error.
  if (any(k_range < 1L))
    stop(sprintf(
      "All values in k_range must be >= 1. Got invalid values: %s",
      paste(sort(unique(k_range[k_range < 1L])), collapse = ", ")
    ))
  cat(sprintf("Running Model Selection across K = %d to %d...\n\n",
              min(k_range), max(k_range)))

  # Resolve a single-type string or a mixed-type named list once, up front, so
  # every K is fit on the same (possibly re-grouped) indicators and descriptor.
  mm          <- .normalize_measurement(measurement, X)
  X           <- mm$indicators
  measurement <- mm$descriptor

  # The same homoscedastic default fit_mixture() resolves, for the same reason
  # (Day 1969; Hathaway 1985). Without it the sweep would compare a different
  # model family than the one a subsequent fit_mixture() call would estimate.
  dots <- list(...)
  if (!"variances_equal" %in% names(dots))
    dots$variances_equal <- is.character(measurement) &&
      length(measurement) == 1L &&
      measurement %in% c("continuous", "continuous_nan",
                         "gaussian_diag", "gaussian_diag_nan")

  results <- list()
  models  <- list()
  for (k in k_range) {
    cat(sprintf("Fitting %d-class model...\n", k))
    fit <- do.call(fit_mixture_internal,
                   c(list(X = X, Y = NULL, n_components = k,
                          measurement = measurement,
                          n_steps = n_steps, n_init = n_init), dots))
    models[[paste0("K", k)]] <- fit
    results[[k]] <- data.frame(
      Classes = k, LL = fit$metrics$ll, Params = fit$metrics$n_params,
      AIC = fit$metrics$aic, BIC = fit$metrics$bic,
      SABIC = fit$metrics$sabic, Entropy = fit$metrics$entropy,
      # No warning is raised anywhere in this loop: it calls the engine
      # directly, and the warning lives in fit_mixture(). The column is how a
      # user learns which K found its maximum only once.
      Unreplicated = .is_unreplicated(fit$metrics)
    )
  }
  fit_table   <- do.call(rbind, results)
  # The rows must be in ascending K for row i to test K against K + 1, which is
  # what k_range gives unless the caller shuffled it.
  fit_table   <- fit_table[order(fit_table$Classes), , drop = FALSE]
  rownames(fit_table) <- NULL
  if (vlmr != "none") {
    cat("Computing the VLMR test...\n")
    fit_table <- .vlmr_augment(fit_table, models, vlmr)
  }
  best_bic_k  <- fit_table$Classes[which.min(fit_table$BIC)]
  cat("\n=== Model Selection Summary ===\n")
  # Rounded column by column: round() on the whole frame fails once one of the
  # columns is logical.
  printed <- fit_table
  num     <- vapply(printed, is.numeric, logical(1))
  printed[num] <- lapply(printed[num], round, 3)
  print(printed)
  flagged <- fit_table$Classes[which(fit_table$Unreplicated)]
  if (length(flagged)) {
    # Every K in one call is fitted at the same `n_init`, so the advice is read
    # off the first flagged model rather than repeated per row.
    m1 <- models[[paste0("K", flagged[1])]]$metrics
    cat(sprintf("\nUnreplicated maximum at K = %s - %s.\n",
                paste(flagged, collapse = ", "),
                .replication_advice(m1$n_requested %||% m1$n_starts)))
  }
  cat(sprintf("\n-> Best model according to BIC: %d classes\n", best_bic_k))
  # Classed purely so that plot() has something to dispatch on. The list is
  # unchanged otherwise, and `result$fit_table` behaves exactly as before.
  out <- list(fit_table = fit_table, models = models, best_k = best_bic_k)
  # The means and standard deviations of the two reference distributions are
  # what tell a suspicious user which distribution produced a p-value, so they
  # are returned; they are not printed, because the table is already wide.
  if (vlmr != "none") {
    out$vlmr <- attr(fit_table, "vlmr_detail")
    attr(out$fit_table, "vlmr_detail") <- NULL
  }
  class(out) <- "mixture_comparison"
  return(out)
}

#' Extract Covariate Odds Ratios from a Fitted Mixture Model
#'
#' @description
#' Extracts the logistic regression coefficients from a covariate structural
#' model and returns them as a matrix of odds ratios, centered on a reference
#' class. Only available when the model was fitted with
#' \code{structural = "covariate"}.
#'
#' @param object A fitted \code{mixture_model} object with a covariate
#'   structural model.
#' @param ref_class Integer. The reference class for centering. All other
#'   class odds ratios are expressed relative to this class. Default is
#'   \code{1}.
#' @param covariate_names Optional character vector of predictor names to
#'   override the column names stored in the model. Default is \code{NULL}.
#' @param exponentiate Logical. When \code{TRUE} (the default) the coefficients
#'   are returned as odds ratios; when \code{FALSE}, as the multinomial-logit
#'   coefficients themselves, relative to \code{ref_class}.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @details
#' The printed summary reports odds ratios because that is the scale these
#' effects are interpreted and published on, and the default here matches it.
#' \code{coef(fit, exponentiate = FALSE)} and
#' \code{\link[=vcov.mixture_model]{vcov}} give the log-scale estimates and
#' their standard errors, for anyone who needs to compare them against another
#' program or pool them across analyses. Both scales are exact —
#' \code{log(coef(fit))} has always recovered the coefficients, since the odds
#' ratios are returned at full double precision; the argument makes that
#' discoverable rather than a trick.
#'
#' @return A K x D numeric matrix, where rows are latent classes and columns
#'   are predictors (including the intercept). With \code{exponentiate = TRUE}
#'   these are odds ratios and the reference class row is all \code{1}; with
#'   \code{exponentiate = FALSE} they are log-odds coefficients and that row is
#'   all \code{0}.
#'
#' @seealso \code{\link[=confint.mixture_model]{confint}} for intervals on the
#'   odds-ratio scale, and \code{\link[=vcov.mixture_model]{vcov}} for the
#'   log-scale standard errors.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' Z <- matrix(rnorm(100), nrow = 100)
#' colnames(Z) <- "age"
#' fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
#'                    structural = "covariate",
#'                    n_steps = 3, correction = "ML", n_init = 5)
#' coef(fit)
#' coef(fit, exponentiate = FALSE)
#'
#' @export
coef.mixture_model <- function(object, ref_class = 1, covariate_names = NULL,
                               exponentiate = TRUE, ...) {
  if (is.null(object$sm) || !inherits(object$sm, "covariate"))
    stop("No covariate model.")
  K     <- object$n_components
  betas <- object$sm$parameters$beta
  if (!is.null(covariate_names))
    colnames(betas) <- c("Intercept", covariate_names)
  betas_ref <- sweep(betas, 2, betas[ref_class, ], "-")
  out       <- if (isTRUE(exponentiate)) exp(betas_ref) else betas_ref
  rownames(out) <- paste("Class", 1:K)
  rownames(out)[ref_class] <- paste("Class", ref_class, "(Ref)")
  return(out)
}
