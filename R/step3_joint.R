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
# Step 1 reaches this objective through the weights alone, so its uncertainty
# is carried by the same first-order term as under ML, with the weights rebuilt
# from the perturbed posteriors in place of the classification table. A lone
# BCH block (a distal outcome without predictors, or class predictors alone) is
# the one-block case of the same objective and takes the same variance.
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
  log_sm <- log_likelihood(if (inherits(sm, "nested")) sm else sm$models[[1L]], Y)
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
# log-likelihood; "robust" is the sandwich around it; "corrected" adds to the
# sandwich the second variance component of Bakk, Oberski & Vermunt (2014, eq.
# 17), J D1 J', the uncertainty carried over from step 1 through the
# classification table (`.nested_step1_term()`). With `Cn = NULL` (BCH) the
# base is the sandwich whatever `se` says, and "corrected" adds the same term
# through the BCH weights; `bch_assign` is then the assigned-class variable
# the weights were built from.
.attach_nested_step3_vcov <- function(model_state, Y, A1, Cn, w, se,
                                      strata = NULL, cluster = NULL,
                                      X = NULL, assignment = "proportional",
                                      bch_assign = NULL) {
  sm   <- model_state$sm
  lone <- !inherits(sm, "nested")
  if (lone) sm <- list(models = stats::setNames(list(sm), class(sm)[1L]))
  if (!all(vapply(sm$models, .nested_packable, logical(1)))) {
    return(model_state)
  }

  pieces <- lapply(sm$models, .nested_block_pack)
  lens   <- lengths(pieces)
  theta  <- unlist(pieces, use.names = FALSE)
  n_par  <- length(theta)
  h      <- 1e-4 * pmax(1, abs(theta))

  # Scores by central differences, then the Hessian by central differences of
  # the summed score -- n_par * 4 evaluations for the scores and the same again
  # for the Hessian, rather than the 2 * n_par^2 of a second difference.
  score_at <- function(par, A = A1, Cv = Cn) {
    S <- matrix(0, nrow(Y), n_par)
    for (j in seq_len(n_par)) {
      e <- replace(numeric(n_par), j, h[j])
      S[, j] <- (.nested_ll_case(sm, Y, A, Cv, par + e, lens) -
                   .nested_ll_case(sm, Y, A, Cv, par - e, lens)) / (2 * h[j])
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

  # A covariate block is fitted with the latent prior, so the joint estimate
  # maximises the step-3 log-likelihood plus that prior: its curvature joins the
  # bread and is taken once out of the meat, as for a lone class regression
  # (the note above .step3_scores(), R/step3_variance.R). The block's free rows
  # are packed in their natural order with the anchor dropped
  # (.nested_block_pack()), which is the helper's order once the anchor is
  # moved last.
  Hp    <- matrix(0, n_par, n_par)
  ends  <- cumsum(lens)
  first <- 1L
  for (j in seq_along(sm$models)) {
    m  <- sm$models[[j]]
    nc <- if (lone) ncol(Y) else sm$columns_per_model[[names(sm$models)[j]]]
    if (inherits(m, "covariate")) {
      B    <- m$parameters$beta
      a    <- .covariate_anchor(B)
      Bord <- sweep(B[c(setdiff(seq_len(nrow(B)), a), a), , drop = FALSE],
                    2, B[a, ], "-")
      Zm   <- .covariate_design(m, Y[, first:(first + nc - 1L), drop = FALSE])
      idx  <- (ends[j] - lens[j] + 1L):ends[j]
      Hp[idx, idx] <- .covariate_prior_derivs(Bord, Zm,
                                              .bayes_alpha(m, "latent"))$hessian
    }
    first <- first + nc
  }
  H <- H + Hp

  B_inv <- pinv(-H)
  if (identical(se, "hessian") && !is.null(Cn)) {
    V      <- B_inv
    method <- "Inverse observed information of the step-3 log-likelihood (joint model)"
  } else {
    meat   <- (if (is.null(strata)) crossprod(S * w)
               else compute_survey_B(S * w, strata, cluster)) - Hp
    V      <- B_inv %*% meat %*% B_inv
    method <- if (is.null(Cn))
      paste("Case-clustered sandwich of the BCH-weighted step-3",
            "log-likelihood (joint model); step-1 uncertainty not propagated")
    else paste("Step-3 sandwich of the joint model; step-1 uncertainty",
               "not propagated")
    if (identical(se, "corrected") && !is.null(X) &&
        (!is.null(Cn) || !is.null(bch_assign))) {
      corr <- .nested_step1_term(model_state, X, A1, w, assignment, B_inv,
                                 function(A, Cv) colSums(w * score_at(theta, A, Cv)),
                                 bch_assign = bch_assign)
      if (is.null(corr)) {
        .warn_step1_unavailable(model_state$mm)
      } else {
        V      <- V + corr$V
        method <- sprintf(paste("Bakk-Oberski-Vermunt corrected (%s%s step 3,",
                                "%s step 1)%s"),
                          if (lone) "" else "joint model; ",
                          if (is.null(Cn)) "BCH case-clustered sandwich"
                          else "robust",
                          corr$method,
                          if (isTRUE(corr$fallback))
                            "; step-1 Hessian was not positive definite" else "")
      }
    }
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
      # A lone block's M-step takes its variances with absolute BCH weights, so
      # they are not where the weighted log-likelihood peaks and their rows of
      # V mean nothing. The means' rows do not depend on them: the mean-variance
      # cross derivatives vanish at the weighted means.
      if (!lone)
        m$parameters$var_ses <- sqrt(pmax(diag(Vj)[K + seq_len(K)], 0))
    } else {
      m$parameters$hessian <- -pinv(Vj)
    }
    sm$models[[j]] <- m
  }

  if (lone) {
    sm <- sm$models[[1L]]
    sm$parameters$V_method <- method
  }
  model_state$sm <- sm
  model_state
}

# The step-1 term of the corrected variance for a joint ML step 3,
# J D1 J' with J = (-H3)^{-1} C and C = d(total step-3 score)/d theta1', as in
# `.step3_covariate_vcov()` but over the joint vector. Moving theta1 moves the
# posteriors and with them the classification table; under proportional
# assignment it moves the assignment weights too. Under modal assignment the
# assignments are held where step 2 put them, since a small change to the
# posteriors must not carry a case across a boundary. `total_score(A, Cv)` is
# the summed step-3 score at the step-3 estimates for assignment `A` and table
# `Cv`. Under BCH (`bch_assign` given) the table is replaced by the BCH weights
# rebuilt from the moved posteriors, with `Cv = NULL`; a modal `bch_assign` is
# held fixed for the same reason. Returns NULL when the measurement model has
# no unconstrained packing.
.nested_step1_term <- function(model_state, X, A1, w, assignment, B_inv,
                               total_score, bch_assign = NULL) {
  th1 <- .step1_pack(model_state)
  if (is.null(th1) || !length(th1)) return(NULL)
  p1 <- length(th1)
  h1 <- .step1_fd_step * pmax(1, abs(th1))

  score_of <- function(v) {
    ms <- .step1_unpack(model_state, v)
    ms$sm <- NULL
    r  <- exp(e_step(ms, X, NULL)$log_resp)
    if (!is.null(bch_assign)) {
      A <- if (identical(assignment, "modal")) bch_assign else r
      return(total_score(.bch_weights(r, A, w), NULL))
    }
    A  <- if (identical(assignment, "modal")) A1 else r
    Cv <- sweep(t(r * w) %*% A, 1, colSums(r * w), "/")
    total_score(A, Cv)
  }
  Cmat <- vapply(seq_len(p1), function(m) {
    a <- th1; a[m] <- a[m] + h1[m]
    b <- th1; b[m] <- b[m] - h1[m]
    (score_of(a) - score_of(b)) / (2 * h1[m])
  }, numeric(nrow(B_inv)))
  dim(Cmat) <- c(nrow(B_inv), p1)

  d1_method <- if (p1 <= .step1_hessian_max) "hessian" else "outer"
  if (d1_method == "outer")
    message(sprintf(paste("Step-1 variance for the corrected standard errors",
                          "uses the outer-product estimator: the measurement",
                          "model has %d parameters (limit %d for the numerical",
                          "Hessian)."), p1, .step1_hessian_max))
  S1 <- .step1_case_scores(model_state, X, th1, w)
  d1 <- .step1_variance(model_state, X, th1, S1, w, d1_method)

  J <- B_inv %*% Cmat
  list(V = J %*% d1$V %*% t(J), method = d1$method, fallback = d1$fallback)
}
