# ==============================================================================
# Standard errors for a joint step-3 model: class predictors and one or more
# distal outcomes held in one nested structural model.
#
# Every block of such a model shares one likelihood -- the class probabilities
# P(k | z) enter every case's contribution alongside each outcome's density --
# so the single-block variance formulas elsewhere do not apply: the covariance
# between a class-regression coefficient and an outcome mean is not zero, and
# dropping it misstates both. The variance here is taken over the whole free
# parameter vector at once, from the step-3 log-likelihood itself,
#
#   l(theta) = sum_i w_i sum_a A1[i, a] log sum_k P(k | z_i) prod_d f_d(y_id | k, z_i) C[k, a],
#
# with the classification-error matrix C and the assigned-class variable A1
# held fixed at their step-1 values (Vermunt, 2010; Bakk, Oberski & Vermunt,
# 2014). The Hessian and the per-case scores are taken by central differences;
# the parameter vector is short (tens of entries) and each evaluation is one
# pass over an n x K matrix, so this is cheap, and a closed form for every
# pairing of blocks would be a great deal of code to keep in step.
#
# Under the BCH correction there is no C matrix: A1 holds the BCH weights d_ik
# and the objective is the weighted log-likelihood
#
#   l(theta) = sum_i w_i sum_k d_ik log [P(k | z_i) prod_d f_d(y_id | k, z_i)],
#
# which is a sum of separate terms and so is maximised by the M-step alone. Its
# variance is always the sandwich with the scores summed per case: the weights
# are not frequencies (many are negative), so the inverse Hessian credits each
# case with a full observation in every class it touches, and the K records of
# one case share its outcome and covariates, so the case is the independent unit.
# ==============================================================================

# Which blocks of a nested structural model this file can pack. A block of any
# other kind has a per-class parameter layout with no cross-class covariance
# slot for the result to be written back to.
.nested_packable <- function(m) {
  inherits(m, c("covariate", "distal_continuous", "distal_continuous_pooled",
                "distal_pooled"))
}

# The anchor row of a covariate block's coefficient matrix: the one held at
# zero. A fresh fit puts it last; a re-sorted one may have moved it, so it is
# found rather than assumed.
.covariate_anchor <- function(beta) {
  zero <- which(rowSums(abs(beta)) == 0)
  if (length(zero)) zero[length(zero)] else nrow(beta)
}

# The free parameters of one block, as a list of (name, value) pieces in a fixed
# order, and the inverse. Variances are packed on their natural scale, which is
# the scale their standard errors are reported on.
.nested_block_pack <- function(m) {
  p <- m$parameters
  if (inherits(m, "covariate")) {
    free <- setdiff(seq_len(nrow(p$beta)), .covariate_anchor(p$beta))
    return(as.vector(t(p$beta[free, , drop = FALSE])))
  }
  if (inherits(m, "distal_continuous_pooled")) {
    v <- if (identical(m$variances, "class_specific")) as.vector(p$covariances)
         else p$covariances[1L, 1L]
    return(c(as.vector(p$beta_pooled), v))
  }
  if (inherits(m, "distal_continuous"))
    return(c(as.vector(p$means), as.vector(p$covariances)))
  # distal_pooled and distal_categorical: row-major, the order the Hessian and
  # .outcome_contrast_parts() index it in.
  as.vector(t(p$beta_pooled))
}

.nested_block_unpack <- function(m, par) {
  p <- m$parameters
  if (inherits(m, "covariate")) {
    free <- setdiff(seq_len(nrow(p$beta)), .covariate_anchor(p$beta))
    p$beta[free, ] <- matrix(par, length(free), ncol(p$beta), byrow = TRUE)
  } else if (inherits(m, "distal_continuous_pooled")) {
    L <- length(p$beta_pooled)
    p$beta_pooled[] <- par[seq_len(L)]
    p$covariances[] <- par[-seq_len(L)]
  } else if (inherits(m, "distal_continuous")) {
    K <- length(p$means)
    p$means[]       <- par[seq_len(K)]
    p$covariances[] <- par[K + seq_len(K)]
  } else {
    p$beta_pooled[] <- matrix(par, nrow(p$beta_pooled), ncol(p$beta_pooled),
                              byrow = TRUE)
  }
  m$parameters <- p
  m
}

# Per-case step-3 log-likelihood at `par`, the quantity the Hessian and the
# scores are taken of. `lens` splits `par` across the blocks in model order.
.nested_ll_case <- function(sm, Y, A1, Cn, par, lens) {
  ends <- cumsum(lens)
  for (j in seq_along(sm$models))
    sm$models[[j]] <- .nested_block_unpack(
      sm$models[[j]], par[(ends[j] - lens[j] + 1L):ends[j]])
  log_sm <- log_likelihood(sm, Y)
  if (is.null(Cn)) return(rowSums(A1 * log_sm))
  lsh    <- .row_max(log_sm)
  Z      <- exp(log_sm - lsh) %*% Cn
  rowSums(A1 * log(pmax(Z, 1e-300))) + lsh
}

# Compute the joint covariance and write each block's share back where that
# block's own readers look for it: V_robust on the covariate block (K*D layout,
# zeros on the anchor), cov_theta/ses/var_ses on a pooled continuous outcome,
# Sigma_mu/ses on a plain continuous one, and the Hessian slot on a categorical
# one, which downstream code inverts. The whole matrix is kept on the nested
# model as `vcov_joint`, for contrasts that cross blocks.
#
# se: "hessian" is the inverse observed information of the step-3
# log-likelihood; "robust" and "corrected" are the sandwich around it. With
# `Cn = NULL` (BCH) the sandwich is reported whatever `se` says. The
# second variance component of Bakk, Oberski & Vermunt (2014), the uncertainty
# carried over from step 1, is not added for a joint model yet, and V_method
# says so.
.attach_nested_step3_vcov <- function(model_state, Y, A1, Cn, w, se,
                                      strata = NULL, cluster = NULL) {
  sm <- model_state$sm
  if (!all(vapply(sm$models, .nested_packable, logical(1)))) {
    return(model_state)
  }

  pieces <- lapply(sm$models, .nested_block_pack)
  lens   <- lengths(pieces)
  theta  <- unlist(pieces, use.names = FALSE)
  n_par  <- length(theta)
  h      <- 1e-4 * pmax(1, abs(theta))

  f  <- function(par) .nested_ll_case(sm, Y, A1, Cn, par, lens)
  fs <- function(par) sum(w * f(par))

  # Scores by central differences, then the Hessian by central differences of
  # the summed score -- n_par * 4 evaluations for the scores and the same again
  # for the Hessian, rather than the 2 * n_par^2 of a second difference.
  score_at <- function(par) {
    S <- matrix(0, nrow(Y), n_par)
    for (j in seq_len(n_par)) {
      e <- replace(numeric(n_par), j, h[j])
      S[, j] <- (f(par + e) - f(par - e)) / (2 * h[j])
    }
    S
  }
  S <- score_at(theta)
  H <- matrix(0, n_par, n_par)
  for (j in seq_len(n_par)) {
    e <- replace(numeric(n_par), j, h[j])
    H[, j] <- (colSums(w * score_at(theta + e)) -
                 colSums(w * score_at(theta - e))) / (2 * h[j])
  }
  H <- (H + t(H)) / 2

  B_inv <- pinv(-H)
  if (identical(se, "hessian") && !is.null(Cn)) {
    V      <- B_inv
    method <- "Inverse observed information of the step-3 log-likelihood (joint model)"
  } else {
    meat   <- if (is.null(strata)) crossprod(S * w)
              else compute_survey_B(S * w, strata, cluster)
    V      <- B_inv %*% meat %*% B_inv
    method <- if (is.null(Cn))
      paste("Case-clustered sandwich of the BCH-weighted step-3",
            "log-likelihood (joint model); step-1 uncertainty not propagated")
    else paste("Step-3 sandwich of the joint model; step-1 uncertainty",
               "not propagated")
  }

  lab <- unlist(lapply(names(sm$models), function(nm)
    paste0(nm, "[", seq_len(lens[[nm]]), "]")), use.names = FALSE)
  dimnames(V) <- list(lab, lab)
  sm$parameters$vcov_joint <- V
  sm$parameters$V_method   <- method

  ends <- cumsum(lens)
  for (j in seq_along(sm$models)) {
    m   <- sm$models[[j]]
    idx <- (ends[j] - lens[j] + 1L):ends[j]
    Vj  <- V[idx, idx, drop = FALSE]
    dimnames(Vj) <- NULL
    if (inherits(m, "covariate")) {
      beta <- m$parameters$beta
      K <- nrow(beta); D <- ncol(beta)
      free <- setdiff(seq_len(K), .covariate_anchor(beta))
      pos  <- as.vector(t(outer((free - 1L) * D, seq_len(D), "+")))
      Vf   <- matrix(0, K * D, K * D)
      Vf[pos, pos] <- Vj
      m$parameters$V_robust <- Vf
      m$parameters$V_method <- method
    } else if (inherits(m, "distal_continuous_pooled")) {
      L <- length(m$parameters$beta_pooled)
      m$parameters$cov_theta <- Vj[seq_len(L), seq_len(L), drop = FALSE]
      m$parameters$ses       <- matrix(sqrt(pmax(diag(Vj)[seq_len(L)], 0)),
                                       nrow = 1)
      m$parameters$var_ses   <- sqrt(pmax(diag(Vj)[-seq_len(L)], 0))
    } else if (inherits(m, "distal_continuous")) {
      K <- length(m$parameters$means)
      m$parameters$Sigma_mu <- Vj[seq_len(K), seq_len(K), drop = FALSE]
      m$parameters$ses      <- matrix(sqrt(pmax(diag(Vj)[seq_len(K)], 0)),
                                      ncol = 1)
      m$parameters$var_ses  <- sqrt(pmax(diag(Vj)[K + seq_len(K)], 0))
    } else {
      m$parameters$hessian <- -pinv(Vj)
    }
    sm$models[[j]] <- m
  }

  model_state$sm <- sm
  model_state
}
