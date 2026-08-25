# ==============================================================================
# mixtureEM Implementation - Utilities and Math Helpers
# ==============================================================================

# Null-coalescing helper. Defined locally so the package carries no dependency
# on rlang or on base R >= 4.4 (which introduced `%||%`).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Log-Sum-Exp trick for numerical stability
logsumexp <- function(x, MARGIN = 1) {
  if (is.null(dim(x))) {
    max_x <- max(x)
    if (max_x == -Inf) return(-Inf)
    return(max_x + log(sum(exp(x - max_x))))
  } else {
    # Row maxima via .row_max() rather than apply(): the row-wise form is the
    # inner loop of the forward-backward recursion in latent transition models
    # and apply() dominates the runtime there.
    max_x <- if (MARGIN == 1) .row_max(x) else apply(x, MARGIN, max)
    max_x[max_x == -Inf] <- 0
    if (MARGIN == 1) {
      res <- max_x + log(rowSums(exp(x - max_x)))
    } else {
      res <- max_x + log(colSums(exp(sweep(x, 2, max_x, "-"))))
    }
    return(res)
  }
}

# Row maxima of a matrix. Same result as apply(x, 1, max), but apply() builds a
# list of rows and calls the closure once per case, which is the single largest
# cost in the E-step of a large fit. max.col() does the search in one pass; the
# tie rule is irrelevant here because only the value is taken, never the index.
.row_max <- function(x) x[cbind(seq_len(nrow(x)), max.col(x, ties.method = "first"))]

# Row indices of cases with no observed value on any measurement indicator.
#
# Such a case carries no information about class membership: under FIML every
# item drops out of its likelihood, so its posterior is exactly the class prior
# and its log-likelihood contribution is exactly zero. It is not free, though —
# it still counts toward the sample size in BIC/SABIC, is assigned modally to
# the largest class, and enters the entropy calculation as a maximally uncertain
# case, which drags relative entropy down. The convention is to delete these
# cases before estimation, and the package follows it.
#
# The pattern is almost always an import artifact (a trailing blank line, a
# mis-parsed delimiter, a merge that did not match) rather than a substantive
# missing-data pattern, which is why the count is reported separately from
# item-level missingness rather than folded into the FIML summary.
.empty_rows <- function(X) {
  if (is.null(X) || nrow(X) == 0L) return(integer(0))
  which(rowSums(!is.na(X)) == 0L)
}

# Subset the cases of a row-aligned argument, whatever shape it arrived in.
# User-supplied covariates and grouping variables may be vectors, factors,
# matrices, or data frames; NULL passes through so callers can apply this
# unconditionally to optional arguments.
.subset_cases <- function(v, keep) {
  if (is.null(v)) return(NULL)
  if (is.null(dim(v))) return(v[keep])
  v[keep, , drop = FALSE]
}

# Render a vector of row indices for a message, keeping it short enough to read.
# Reporting which rows were dropped is what lets a user trace the problem back
# to the file they imported, so the first few are always named.
.abbreviate_indices <- function(idx, max_show = 6L) {
  if (length(idx) <= max_show) return(paste(idx, collapse = ", "))
  paste0(paste(idx[seq_len(max_show)], collapse = ", "),
         ", ... (", length(idx) - max_show, " more)")
}

# Shorten display labels that would not fit a printed table column. Labels at
# or under `width` characters pass through unchanged; longer ones are
# abbreviated to `width` (with a strtrim fallback for non-ASCII text, which
# abbreviate() refuses) and disambiguated when the abbreviation collides.
# The mapping from shortened to full label is attached as attr(, "legend") so
# the caller can print a key; .cat_label_legend() does exactly that.
# Data frames returned to the user always carry the full names — shortening
# is a display concern only.
.shorten_labels <- function(x, width = 28L) {
  x <- as.character(x)
  too_long <- !is.na(x) & nchar(x) > width
  if (!any(too_long)) return(x)

  short <- vapply(x[too_long], function(s) {
    if (grepl("[^\x01-\x7F]", s)) paste0(strtrim(s, width - 1L), "~")
    else abbreviate(s, minlength = width, strict = TRUE, dot = FALSE)
  }, character(1), USE.NAMES = FALSE)

  # Uniqueness against both the unchanged labels and each other; make.unique
  # can push a label a couple of characters past `width`, which beats showing
  # two identical rows for different variables.
  short <- make.unique(c(x[!too_long], short), sep = "~")
  short <- short[(sum(!too_long) + 1L):length(x)]

  out <- x
  out[too_long] <- short
  attr(out, "legend") <- stats::setNames(x[too_long], short)
  out
}

# Width of the label column for a printed table: wide enough for the longest
# (already shortened) label, never narrower than `min`.
.label_width <- function(x, min = 13L) {
  if (!length(x)) return(as.integer(min))
  max(as.integer(min), max(nchar(x), 0L, na.rm = TRUE))
}

# Print the key for labels shortened by .shorten_labels(), one per line.
# Prints nothing when no label was shortened.
.cat_label_legend <- function(labels, indent = "  ") {
  legend <- attr(labels, "legend")
  if (is.null(legend) || !length(legend)) return(invisible(NULL))
  cat(sprintf("%sAbbreviated names:\n", indent))
  for (i in seq_along(legend))
    cat(sprintf("%s  %s = %s\n", indent, names(legend)[i], legend[[i]]))
  invisible(NULL)
}

# Line up a block of fitted item parameters with the columns of the stored
# indicator matrix, so that anything reported per item can be put beside the
# data it came from.
#
# By name first, which is the only rule that can be trusted: a mixed measurement
# model fits its blocks in an order that need not be the data's, and a sorted or
# subset fit may have moved them again. Positionally only when there are no
# names to go on *and* the block accounts for every column of the data, which is
# the single-block unnamed-matrix case where the fitting order is the data
# order by construction. NULL when neither rule applies - the caller must then
# do without rather than guess, since a wrong alignment reports one item's
# number under another item's name.
.match_indicator_columns <- function(cols, data) {
  if (is.null(cols) || is.null(data) || !length(cols)) return(NULL)
  matched <- match(cols, colnames(data))
  if (anyNA(matched) && length(cols) == ncol(data) && is.null(colnames(data)))
    matched <- seq_along(cols)
  if (anyNA(matched)) return(NULL)
  matched
}

# Resolve case weights and the effective sample size used by BIC and SABIC.
#
# Two kinds of weight mean different things:
#
#   "sampling"  -- probability or design weights, saying how much of the
#                  population each case stands for. Only their relative sizes
#                  carry information, so they are rescaled to sum to the number
#                  of observed cases. This keeps the effective sample size equal
#                  to the number of cases actually collected, which is what
#                  stops BIC from being inflated by an arbitrary weight scale.
#
#   "frequency" -- counts of identical cases, as in a response-pattern table
#                  where one row stands for many respondents. Here the total
#                  really is the sample size, so the weights are used unchanged
#                  and the effective sample size is their sum.
.resolve_weights <- function(weights, n_samples, weight_type = "sampling",
                             label = "weights") {
  weight_type <- match.arg(weight_type, c("sampling", "frequency"))

  if (is.null(weights))
    return(list(weights = rep(1, n_samples), n_eff = n_samples,
                type = weight_type))

  weights <- as.numeric(weights)
  if (length(weights) != n_samples)
    stop(sprintf("`%s` must have one entry per case (expected %d, got %d).",
                 label, n_samples, length(weights)), call. = FALSE)
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights < 0))
    stop(sprintf("`%s` must be finite and non-negative.", label), call. = FALSE)
  if (sum(weights) <= 0)
    stop(sprintf("`%s` must include at least one positive value.", label),
         call. = FALSE)

  if (weight_type == "frequency")
    return(list(weights = weights, n_eff = sum(weights), type = "frequency"))

  list(weights = weights / sum(weights) * n_samples,
       n_eff = n_samples, type = "sampling")
}

# Do the supplied weights look like frequency counts rather than sampling
# weights? Whole numbers that add up to far more than the number of rows are the
# signature of a response-pattern table. Used only to point out a likely
# mis-specification, never to override what was asked for.
.looks_like_frequencies <- function(weights, n_samples) {
  if (is.null(weights)) return(FALSE)
  w <- as.numeric(weights)
  if (anyNA(w) || any(!is.finite(w))) return(FALSE)
  all(abs(w - round(w)) < 1e-8) && all(w >= 1) && sum(w) > 1.5 * n_samples
}

# Row-wise softmax of a matrix of logits, stabilised by the row maximum.
# max.col() rather than apply() because this sits inside the inner loop of the
# multinomial-logit fitter, which an EM algorithm calls once per iteration.
softmax_rows <- function(logits, clamp = 50) {
  logits <- pmax(pmin(logits, clamp), -clamp)
  max_l  <- logits[cbind(seq_len(nrow(logits)),
                         max.col(logits, ties.method = "first"))]
  e <- exp(logits - max_l)
  e / rowSums(e)
}

# Input validation: Check if positive
check_positive <- function(...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!is.numeric(val) || any(val <= 0, na.rm = TRUE))
      stop(sprintf("Expected %s > 0, got %s", name, paste(val, collapse=" ")), call. = FALSE)
  }
}

# Input validation: Check if non-negative
check_nonneg <- function(...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!is.numeric(val) || any(val < 0, na.rm = TRUE))
      stop(sprintf("Expected %s >= 0, got %s", name, paste(val, collapse=" ")), call. = FALSE)
  }
}

# Input validation: Check if value is in allowed choices
check_in <- function(choices, ...) {
  args <- list(...)
  for (name in names(args)) {
    val <- args[[name]]
    if (!all(val %in% choices))
      stop(sprintf("%s value '%s' not recognized. Choose from: %s",
                   name, paste(val, collapse=" "), paste(choices, collapse=", ")),
           call. = FALSE)
  }
}

# Modal (Hard) Assignment of Probabilities
modal <- function(resp, clip = FALSE) {
  max_idx <- max.col(resp, ties.method = "random")
  modal_resp <- matrix(0, nrow = nrow(resp), ncol = ncol(resp))
  modal_resp[cbind(1:nrow(resp), max_idx)] <- 1
  if (clip) modal_resp <- pmax(pmin(modal_resp, 1 - 1e-15), 1e-15)
  return(modal_resp)
}

# Multiple categorical one-hot encoding (0-indexed input)
max_one_hot <- function(mat, max_n_outcomes = NULL, total_outcomes = NULL) {
  n_samples <- nrow(mat)
  n_features <- ncol(mat)

  if (is.null(max_n_outcomes) || is.null(total_outcomes)) {
    outcomes <- apply(mat, 2, max, na.rm = TRUE) + 1
    total_outcomes <- sum(outcomes)
    max_n_outcomes <- max(outcomes)
  }

  one_hot <- matrix(0, nrow = n_samples, ncol = n_features * max_n_outcomes)

  for (c in 1:n_features) {
    integer_codes <- mat[, c]
    not_observed <- is.na(integer_codes)
    integer_codes[not_observed] <- 0
    col_indices <- integer_codes + (c - 1) * max_n_outcomes + 1
    one_hot[cbind(1:n_samples, col_indices)] <- 1
    if (any(not_observed)) {
      start_col <- (c - 1) * max_n_outcomes + 1
      end_col   <- c * max_n_outcomes
      one_hot[not_observed, start_col:end_col] <- NA
    }
  }

  return(list(one_hot = one_hot,
              max_n_outcomes = max_n_outcomes,
              total_outcomes = total_outcomes))
}

# Moore-Penrose pseudoinverse (pure base R)
pinv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X)
  d <- s$d
  d_inv <- ifelse(d > tol * d[1], 1/d, 0)
  return(s$v %*% diag(d_inv, nrow = length(d)) %*% t(s$u))
}

# Design-adjusted "meat" matrix for the robust sandwich estimator under a
# complex survey design (Taylor-series linearization).
#
# Given an N x P matrix of individual-level score (gradient) vectors, the
# scores are aggregated to the primary sampling unit (PSU) level within each
# stratum and combined as
#
#   B = sum_o (C_o / (C_o - 1)) sum_c (g_oc - g_bar_o)(g_oc - g_bar_o)^T
#
# where o indexes strata, c indexes the C_o PSUs in stratum o, g_oc is the
# summed score over all cases in PSU c, and g_bar_o is the mean PSU score in
# the stratum. The finite population correction is omitted, which is the
# standard choice for large-scale surveys where it approaches 1.
#
# PSUs are identified within strata: an identical cluster label appearing in
# two different strata denotes two distinct PSUs.
#
# Singleton strata (C_o = 1) make the C_o/(C_o - 1) multiplier undefined.
# They are handled with the "adjust" / centered-at-grand-mean convention:
# the single PSU is centered on the overall mean of all PSU scores and the
# multiplier defaults to 1, which contributes a conservative variance term
# rather than dropping the stratum.
compute_survey_B <- function(score_mat, strata, cluster) {
  score_mat <- as.matrix(score_mat)
  P <- ncol(score_mat)
  B <- matrix(0, nrow = P, ncol = P)

  # Combined key so PSUs are unique within (and only within) their stratum.
  psu_key <- paste(strata, cluster, sep = "\r")

  # Grand mean of PSU-aggregated scores, used to center singleton strata.
  unique_psu <- unique(psu_key)
  g_all_psu  <- matrix(0, nrow = length(unique_psu), ncol = P)
  for (i in seq_along(unique_psu)) {
    rows <- which(psu_key == unique_psu[i])
    g_all_psu[i, ] <- colSums(score_mat[rows, , drop = FALSE])
  }
  g_bar_grand <- colMeans(g_all_psu)

  for (o in unique(strata)) {
    in_stratum     <- which(strata == o)
    psu_in_stratum <- unique(psu_key[in_stratum])
    C_o            <- length(psu_in_stratum)

    # Aggregate individual scores to the PSU level within this stratum.
    g_oc <- matrix(0, nrow = C_o, ncol = P)
    for (c_idx in seq_len(C_o)) {
      rows         <- in_stratum[psu_key[in_stratum] == psu_in_stratum[c_idx]]
      g_oc[c_idx, ] <- colSums(score_mat[rows, , drop = FALSE])
    }

    if (C_o == 1L) {
      # Singleton stratum: center on the grand mean, multiplier = 1.
      centered  <- sweep(g_oc, 2, g_bar_grand, "-")
      stratum_B <- crossprod(centered)
    } else {
      # Standard linearization: center on the stratum mean.
      g_bar_o   <- colMeans(g_oc)
      centered  <- sweep(g_oc, 2, g_bar_o, "-")
      stratum_B <- (C_o / (C_o - 1)) * crossprod(centered)
    }

    B <- B + stratum_B
  }

  return(B)
}

# Relative Entropy (normalised by log(K), scales 0-1)
relative_entropy <- function(absolute_entropy, n_samples, n_classes) {
  if (n_classes <= 1) return(1)
  rel_ent <- 1 - (absolute_entropy / (n_samples * log(n_classes)))
  return(max(0, min(1, rel_ent)))
}

# Completion of incomplete covariates for structural models under the
# endogenous-constrained-x approach (Sterba, 2014, Multivariate Behavioral
# Research, 49, 614-632; see also Depaoli, Jia & Visser, 2025).
#
# Rather than listwise deleting cases with missing covariates or filling them
# with an unconditional column mean, the covariates are treated as endogenous
# with a single Gaussian marginal distribution whose parameters are held
# common across latent classes. Each missing block of a row is filled with its
# conditional expectation given the observed block under that shared Gaussian,
#
#   E[x_mis | x_obs] = mu_mis + Sigma_{mis,obs} Sigma_{obs,obs}^{-1} (x_obs - mu_obs),
#
# which is the best linear predictor and uses the inter-covariate associations
# (the off-diagonal of Sigma) instead of ignoring them. Holding the marginal
# class-invariant is what lets the joint likelihood recover the conditional
# model: the marginal density cancels from the class posteriors (Sterba, 2014,
# eq. 13), so completion never lets covariate shape define the class structure.
#
# On complete data this function is an exact no-op, preserving the equivalence
# between the conditional (exogenous-x) and joint (endogenous-constrained-x)
# specifications that holds when covariates are fully observed.
complete_covariates <- function(Z) {
  Z <- as.matrix(Z)
  if (ncol(Z) == 0L || !anyNA(Z)) return(Z)

  mu <- colMeans(Z, na.rm = TRUE)

  # A column with no observed values cannot inform any conditional mean. Centre
  # it at zero and warn, matching the previous all-NA safeguard.
  all_na <- is.nan(mu)
  if (any(all_na)) {
    warning(sprintf(
      paste0("complete_covariates: column(s) %s are entirely NA. ",
             "Centring at 0. Check your data for completely missing covariates."),
      paste(which(all_na), collapse = ", ")))
    mu[all_na] <- 0
  }

  p <- ncol(Z)
  if (p == 1L) {
    # A single covariate has nothing to condition on, so the conditional mean
    # under the shared Gaussian reduces to the marginal mean.
    na_idx <- is.na(Z[, 1L])
    Z[na_idx, 1L] <- mu[1L]
    return(Z)
  }

  # Class-invariant covariance of the covariates. Pairwise-complete estimation
  # uses every observed cell; a small ridge keeps the observed blocks solvable.
  Sigma <- stats::cov(Z, use = "pairwise.complete.obs")
  Sigma[is.na(Sigma)] <- 0
  diag(Sigma) <- diag(Sigma) + 1e-6

  miss_rows <- which(rowSums(is.na(Z)) > 0L)
  for (i in miss_rows) {
    mis <- which(is.na(Z[i, ]))
    obs <- which(!is.na(Z[i, ]))
    if (length(obs) == 0L) {
      Z[i, mis] <- mu[mis]                 # whole row missing: marginal mean
      next
    }
    S_oo_inv  <- pinv(Sigma[obs, obs, drop = FALSE])
    S_mo      <- Sigma[mis, obs, drop = FALSE]
    Z[i, mis] <- mu[mis] +
      as.vector(S_mo %*% S_oo_inv %*% (Z[i, obs] - mu[obs]))
  }
  Z
}

# Complete only the covariate columns carried by a structural sub-model, under
# the endogenous-constrained-x scheme. For distal-outcome models the first
# column of Y is the endogenous outcome and is left untouched (its own
# likelihood masks missing values via FIML); the remaining columns are
# covariates/moderators. For a class-prediction (covariate) model every column
# is a covariate.
complete_structural_covariates <- function(sm, Y) {
  if (is.null(Y) || is.null(sm)) return(Y)
  Y <- as.matrix(Y)
  distal_types <- c("distal_continuous", "distal_continuous_regression",
                    "distal_continuous_pooled", "distal_pooled",
                    "distal_regression", "distal_categorical")
  if (inherits(sm, distal_types)) {
    if (ncol(Y) > 1L)
      Y[, -1L] <- complete_covariates(Y[, -1L, drop = FALSE])
  } else {
    Y <- complete_covariates(Y)
  }
  Y
}

# ==============================================================================
# Shared helpers for categorical distal outcome models
# (used by both distal_regression.R and distal_pooled.R)
# ==============================================================================

# One-hot encode a 1-indexed integer vector Y into an (n x max_val) matrix.
# Missing values produce all-zero rows.
distal_one_hot <- function(Y, max_val) {
  n   <- length(Y)
  out <- matrix(0, nrow = n, ncol = max_val)
  valid <- !is.na(Y)
  if (any(valid)) out[cbind(which(valid), Y[valid])] <- 1
  return(out)
}

# Multinomial softmax forward pass.
# The first category is the reference (fixed logit = 0).
# @param Z        N x D design matrix
# @param beta_matrix (M-1) x D coefficient matrix
distal_forward <- function(Z, beta_matrix) {
  logits_active <- Z %*% t(beta_matrix)
  logits_full   <- cbind(0, logits_active)
  max_logits    <- .row_max(logits_full)
  exp_logits    <- exp(sweep(logits_full, 1, max_logits, "-"))
  return(sweep(exp_logits, 1, rowSums(exp_logits), "/"))
}

# ==============================================================================
# prepare_covariates()
#
# Converts a user-supplied Y (data.frame, matrix, or vector) into a plain
# numeric matrix suitable for structural models, with two behaviours:
#
#   * Numeric / integer columns are passed through unchanged.
#   * Factor / character columns are dummy-coded (reference = first level):
#       - 2-level factor  → 1 binary column named after the variable
#       - k-level factor  → k-1 columns named "var.level2", "var.level3", ...
#
# Column names are always preserved / generated, so downstream summary
# functions can display real variable names instead of Z1, Z2, etc.
#
# The result also carries a "covariate_terms" attribute: one entry per column
# naming the *original* variable it came from, so that the k-1 dummies of a
# k-level factor can be recognised as one term. Column names alone cannot do
# this — recovering the grouping from "Marital.Single" by pattern matching is
# ambiguous the moment two variables share a prefix (Age and Age_Decades) — and
# the grouping is what an omnibus test of a multi-category covariate needs.
# ==============================================================================
prepare_covariates <- function(Y, min_level_n = 5L) {
  # Plain numeric matrix — nothing to do. An already-prepared matrix passes
  # through here on its way into fit_mixture_internal(), so returning it
  # untouched is also what preserves the terms attribute set below.
  if (is.matrix(Y) && is.numeric(Y)) return(Y)

  # Bare numeric vector — wrap in single-column matrix
  if (is.numeric(Y) && is.null(dim(Y)))
    return(matrix(Y, ncol = 1L,
                  dimnames = list(NULL, deparse(substitute(Y)))))

  df <- as.data.frame(Y)
  out_cols <- vector("list", ncol(df))
  out_terms <- vector("list", ncol(df))
  idx <- 0L

  for (nm in names(df)) {
    col <- df[[nm]]

    # ── numeric / integer: pass through ──────────────────────────────────────
    if (is.numeric(col) || is.integer(col)) {
      idx <- idx + 1L
      m   <- matrix(as.numeric(col), ncol = 1L,
                    dimnames = list(NULL, nm))
      out_cols[[idx]]  <- m
      out_terms[[idx]] <- nm
      next
    }

    # ── factor / character: dummy-code ───────────────────────────────────────
    if (!is.factor(col)) col <- factor(col)

    # Levels with no observations get no dummy. They typically appear when the
    # data were subset (e.g. one country from a pooled file) without
    # droplevels(); keeping them would print a degenerate OR = 1.00 [1, 1]
    # row that reads like a result but is an artifact of the empty column.
    n_unused <- sum(table(col) == 0L)
    if (n_unused > 0L) {
      col <- droplevels(col)
      message(sprintf(
        "prepare_covariates: dropped %d unused level%s of '%s'.",
        n_unused, if (n_unused == 1L) "" else "s", nm))
    }
    lvls <- levels(col)

    if (length(lvls) < 2L) {
      warning(sprintf(
        "prepare_covariates: variable '%s' has fewer than 2 levels and will be dropped.",
        nm))
      next
    }

    # A level observed a handful of times still gets its dummy — dropping
    # observed data silently would be worse — but its coefficient rests on
    # almost no information, so say so before the summary prints an enormous
    # confidence interval without explanation.
    counts <- table(col)
    small  <- counts[counts < min_level_n]
    if (length(small))
      warning(sprintf(
        paste0("Covariate '%s': level%s %s ha%s fewer than %d observations; ",
               "the corresponding coefficients will be unstable. Consider ",
               "merging sparse categories."),
        nm, if (length(small) == 1L) "" else "s",
        paste0("'", names(small), "' (n = ", small, ")", collapse = ", "),
        if (length(small) == 1L) "s" else "ve", min_level_n),
        call. = FALSE)

    other <- lvls[-1L]   # reference = first level

    if (length(other) == 1L) {
      # Binary: single column named "var.Level" so the active level is clear
      idx <- idx + 1L
      out_cols[[idx]] <- matrix(as.integer(col == other),
                                ncol = 1L,
                                dimnames = list(NULL, paste0(nm, ".", other)))
      out_terms[[idx]] <- nm
    } else {
      # Multicategorical: k-1 columns named "var.levelX"
      idx <- idx + 1L
      m   <- matrix(0L, nrow = length(col), ncol = length(other),
                    dimnames = list(NULL, paste0(nm, ".", other)))
      for (j in seq_along(other))
        m[, j] <- as.integer(col == other[j])
      out_cols[[idx]]  <- m
      out_terms[[idx]] <- rep(nm, length(other))
    }
  }

  out <- do.call(cbind, out_cols[seq_len(idx)])
  attr(out, "covariate_terms") <- unlist(out_terms[seq_len(idx)], use.names = FALSE)
  out
}

# ==============================================================================
# Term grouping for a prepared covariate matrix
# ==============================================================================
#
# `.covariate_terms()` answers "which original variable did each column come
# from?". It reads the attribute prepare_covariates() sets, and falls back to
# the column names when there is none — which is the right answer for a plain
# numeric matrix the user supplied directly, where every column is its own term.
#
# `.cbind_covariates()` is cbind() that carries the grouping across, since
# cbind() drops attributes; the covariate matrix is assembled in pieces when an
# outcome column or a grouping variable's design is bound alongside it.
.covariate_terms <- function(M) {
  if (is.null(M)) return(NULL)
  tm <- attr(M, "covariate_terms")
  if (!is.null(tm) && length(tm) == ncol(M)) return(as.character(tm))
  nms <- colnames(M)
  if (!is.null(nms)) return(as.character(nms))
  paste0("V", seq_len(ncol(M)))
}

.cbind_covariates <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(NULL)
  terms <- unlist(lapply(parts, .covariate_terms), use.names = FALSE)
  out   <- do.call(cbind, parts)
  attr(out, "covariate_terms") <- terms
  out
}
