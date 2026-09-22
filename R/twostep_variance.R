# ==============================================================================
# Standard errors for the two-step estimator
# ==============================================================================
#
# A two-step fit (R/stepwise.R, `n_steps == 2`) holds the measurement model at
# its step-1 estimate theta1 and maximises the joint log-likelihood
#
#     L(theta1, theta2) = sum_i w_i log sum_k P(k | z_i; theta2) P(x_i | k; theta1)
#
# over the class-prediction coefficients theta2 alone. Treating theta1 as known
# understates the variance of theta2, because theta1 is itself an estimate
# whose sampling error propagates into the second step. The pseudo-maximum-
# likelihood variance (Gong & Samaniego, 1981), written for this estimator by
# Bakk & Kuha (2018, eq. 5), is
#
#     V = V2 + V1 ,   V2 = I22^{-1} ,   V1 = I22^{-1} I12' S11 I12 I22^{-1}
#
# with I the observed information of the joint log-likelihood at (theta1,
# theta2), partitioned by parameter block, and S11 the sampling variance of the
# step-1 estimate. V2 is what an analysis that takes theta1 as given reports;
# V1 is the part it omits, and it is not small when the classes are poorly
# separated (Bakk & Kuha, Table 2).
#
# I22 and I12 are differenced numerically from the joint log-likelihood. The
# measurement block enters that likelihood only through log P(x_i | k), which
# is fixed while theta2 moves, and the structural block only through
# log P(k | z_i), which is fixed while theta1 moves, so the two blocks are
# differenced separately rather than through one Hessian over the concatenated
# vector: I22 re-evaluates the structural densities alone, and I12 combines
# one pair of perturbed measurement densities per theta1 coordinate with one
# pair of perturbed structural densities per theta2 coordinate. S11 comes from
# .step1_variance() (R/step3_variance.R), the same estimator the three-step
# correction uses, restricted to the measurement-parameter block: with class
# predictors, theta1 is the item parameters only -- the step-1 class sizes
# are discarded and re-estimated as the intercepts of the predictor model
# (Bakk & Kuha, sec. 2.3) -- and the variance of that block is the
# corresponding block of the inverse step-1 information, not the inverse of
# the block.
#
# References
#   Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#     classes and external variables. Psychometrika, 83(4), 871-892.
#   Gong, G., & Samaniego, F. J. (1981). Pseudo maximum likelihood estimation:
#     Theory and applications. The Annals of Statistics, 9(4), 861-869.

# The blocks of the joint information matrix a two-step covariate fit needs,
# by central differences on the relative step every step-one derivative uses
# (`.step1_fd_step`). Returns H22 (p2 x p2) and H12 (p1 x p2), both as second
# derivatives of the log-likelihood (negate for information), and the n x p2
# matrix of case-level scores in theta2. H12 is NULL when the measurement
# model has no unconstrained packing.
.twostep_information <- function(model_state, X, Y, w) {
  mm <- model_state$mm
  sm <- model_state$sm
  th2 <- .step1_pack_sm(sm)
  p2  <- length(th2)
  h2  <- .step1_fd_step * pmax(1, abs(th2))

  A0   <- log_likelihood(mm, X)
  B_at <- function(v) log_likelihood(.step1_unpack_sm(sm, v), Y)
  ll   <- function(A, B) sum(w * logsumexp(A + B, MARGIN = 1))

  # One pair of perturbed structural densities per theta2 coordinate; the
  # scores and the cross block both read from these.
  Bp <- lapply(seq_len(p2), function(b) { v <- th2; v[b] <- v[b] + h2[b]; B_at(v) })
  Bm <- lapply(seq_len(p2), function(b) { v <- th2; v[b] <- v[b] - h2[b]; B_at(v) })

  S2 <- vapply(seq_len(p2), function(b)
    w * (logsumexp(A0 + Bp[[b]], MARGIN = 1) -
           logsumexp(A0 + Bm[[b]], MARGIN = 1)) / (2 * h2[b]),
    numeric(nrow(X)))
  dim(S2) <- c(nrow(X), p2)

  # H22: the four-point second difference of .step1_fd_hessian(), with the
  # measurement densities held fixed.
  H22 <- matrix(0, p2, p2)
  for (a in seq_len(p2)) {
    for (b in a:p2) {
      ea <- numeric(p2); ea[a] <- h2[a]
      eb <- numeric(p2); eb[b] <- h2[b]
      H22[a, b] <- H22[b, a] <-
        (ll(A0, B_at(th2 + ea + eb)) - ll(A0, B_at(th2 + ea - eb)) -
           ll(A0, B_at(th2 - ea + eb)) + ll(A0, B_at(th2 - ea - eb))) /
        (4 * h2[a] * h2[b])
    }
  }

  th1 <- .step1_pack_mm(mm)
  if (is.null(th1) || !length(th1))
    return(list(H22 = H22, H12 = NULL, S2 = S2))

  p1  <- length(th1)
  h1  <- .step1_fd_step * pmax(1, abs(th1))
  H12 <- matrix(0, p1, p2)
  for (a in seq_len(p1)) {
    va <- th1; va[a] <- va[a] + h1[a]
    vb <- th1; vb[a] <- vb[a] - h1[a]
    Ap <- log_likelihood(.step1_unpack_mm(mm, va), X)
    Am <- log_likelihood(.step1_unpack_mm(mm, vb), X)
    for (b in seq_len(p2)) {
      H12[a, b] <- (ll(Ap, Bp[[b]]) - ll(Ap, Bm[[b]]) -
                      ll(Am, Bp[[b]]) + ll(Am, Bm[[b]])) / (4 * h1[a] * h2[b])
    }
  }
  list(H22 = H22, H12 = H12, S2 = S2)
}

# Variance-covariance matrix of a two-step covariate model.
#
# @param model_state Fitted two-step model: `mm` and `weights` at their step-1
#   values, `sm` the covariate model at the joint maximum.
# @param X Indicator matrix the measurement model was fitted to.
# @param Y Covariate data the structural model was fitted to.
# @param w Case weights.
# @param se One of "corrected", "robust", "hessian".
#
# @return A list with `V` (the K*D square covariance matrix, the anchored class
#   padded with zeros) and `method`, a label naming the estimator used.
.twostep_covariate_vcov <- function(model_state, X, Y, w, se = "corrected") {
  K  <- model_state$n_components
  sm <- model_state$sm
  D  <- ncol(sm$parameters$beta)

  # The anchor is row K here as in .step3_covariate_vcov(): every caller runs
  # before sort_model_classes(), which re-anchors beta and this matrix together.
  pad <- function(V) {
    out <- matrix(0, K * D, K * D)
    if ((K - 1L) * D > 0L)
      out[seq_len((K - 1L) * D), seq_len((K - 1L) * D)] <- V
    out
  }
  if (K < 2L || D < 1L)
    return(list(V = matrix(0, K * D, K * D), method = "None (no free parameters)"))

  info <- .twostep_information(model_state, X, Y, w)
  V2   <- pinv(-info$H22)
  V2   <- (V2 + t(V2)) / 2

  if (identical(se, "hessian"))
    return(list(V = pad(V2), method = "Two-step observed information (step 2 only)"))

  if (identical(se, "robust")) {
    strata <- if (isTRUE(sm$has_survey_design) && length(sm$strata) == nrow(X))
      sm$strata else rep(1L, nrow(X))
    cluster <- if (isTRUE(sm$has_survey_design) && length(sm$cluster) == nrow(X))
      sm$cluster else seq_len(nrow(X))
    D2 <- V2 %*% compute_survey_B(info$S2, strata, cluster) %*% V2
    return(list(V = pad((D2 + t(D2)) / 2),
                method = sprintf("Two-step sandwich (%s)",
                                 if (isTRUE(sm$has_survey_design))
                                   "survey-linearized" else "robust")))
  }

  # --- The step-one term ------------------------------------------------------
  th1 <- .step1_pack(model_state)
  if (is.null(info$H12) || is.null(th1) || !length(th1))
    return(list(V = pad(V2),
                method = paste("Two-step observed information (step 2 only);",
                               "step-1 term unavailable for this measurement model")))

  p1 <- length(th1)
  d1_method <- if (p1 <= .step1_hessian_max) "hessian" else "outer"
  if (d1_method == "outer")
    message(sprintf(paste("Step-1 variance for the two-step standard errors",
                          "uses the outer-product estimator: the measurement",
                          "model has %d parameters (limit %d for the numerical",
                          "Hessian)."), p1, .step1_hessian_max))
  S1  <- .step1_case_scores(model_state, X, th1, w)
  d1  <- .step1_variance(model_state, X, th1, S1, w, d1_method)
  nm  <- nrow(info$H12)
  S11 <- d1$V[seq_len(nm), seq_len(nm), drop = FALSE]

  V1 <- V2 %*% t(info$H12) %*% S11 %*% info$H12 %*% V2
  V  <- V2 + V1
  V  <- (V + t(V)) / 2

  list(V = pad(V),
       method = sprintf("Two-step pseudo-ML (Bakk and Kuha, 2018; %s step 1)%s",
                        d1$method,
                        if (isTRUE(d1$fallback))
                          "; step-1 Hessian was not positive definite" else ""))
}

# Compute and attach the two-step variance to a covariate sub-model, the way
# .attach_step3_covariate_vcov() does for the third step. A structural model
# that is not a class-prediction one (a distal outcome) has no unconstrained
# packing here and keeps the Q-function Hessian its readers already label.
.attach_twostep_covariate_vcov <- function(model_state, X, Y, w, se = "corrected") {
  if (!inherits(model_state$sm, "covariate")) return(model_state)
  if (identical(se, "none")) return(model_state)

  res <- tryCatch(
    .twostep_covariate_vcov(model_state, X, Y, w, se = se),
    error = function(e) {
      warning(sprintf(
        "Two-step covariate standard errors fell back to the Q-function Hessian: %s",
        conditionMessage(e)), call. = FALSE)
      NULL
    })
  if (is.null(res)) return(model_state)

  model_state$sm$parameters$V_robust <- res$V
  model_state$sm$parameters$V_method <- res$method
  model_state
}
