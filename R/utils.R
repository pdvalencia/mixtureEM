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
    max_x <- apply(x, MARGIN, max)
    max_x[max_x == -Inf] <- 0
    if (MARGIN == 1) {
      res <- max_x + log(rowSums(exp(x - max_x)))
    } else {
      res <- max_x + log(colSums(exp(sweep(x, 2, max_x, "-"))))
    }
    return(res)
  }
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
  max_logits    <- apply(logits_full, 1, max)
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
# ==============================================================================
prepare_covariates <- function(Y) {
  # Plain numeric matrix — nothing to do
  if (is.matrix(Y) && is.numeric(Y)) return(Y)

  # Bare numeric vector — wrap in single-column matrix
  if (is.numeric(Y) && is.null(dim(Y)))
    return(matrix(Y, ncol = 1L,
                  dimnames = list(NULL, deparse(substitute(Y)))))

  df <- as.data.frame(Y)
  out_cols <- vector("list", ncol(df))
  idx <- 0L

  for (nm in names(df)) {
    col <- df[[nm]]

    # ── numeric / integer: pass through ──────────────────────────────────────
    if (is.numeric(col) || is.integer(col)) {
      idx <- idx + 1L
      m   <- matrix(as.numeric(col), ncol = 1L,
                    dimnames = list(NULL, nm))
      out_cols[[idx]] <- m
      next
    }

    # ── factor / character: dummy-code ───────────────────────────────────────
    if (!is.factor(col)) col <- factor(col)
    lvls <- levels(col)

    if (length(lvls) < 2L) {
      warning(sprintf(
        "prepare_covariates: variable '%s' has fewer than 2 levels and will be dropped.",
        nm))
      next
    }

    other <- lvls[-1L]   # reference = first level

    if (length(other) == 1L) {
      # Binary: single column named "var.Level" so the active level is clear
      idx <- idx + 1L
      out_cols[[idx]] <- matrix(as.integer(col == other),
                                ncol = 1L,
                                dimnames = list(NULL, paste0(nm, ".", other)))
    } else {
      # Multicategorical: k-1 columns named "var.levelX"
      idx <- idx + 1L
      m   <- matrix(0L, nrow = length(col), ncol = length(other),
                    dimnames = list(NULL, paste0(nm, ".", other)))
      for (j in seq_along(other))
        m[, j] <- as.integer(col == other[j])
      out_cols[[idx]] <- m
    }
  }

  do.call(cbind, out_cols[seq_len(idx)])
}
