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
# the block. With a distal outcome the class sizes are not re-estimated: they
# enter P(k) in place of P(k | z_i), so theta1 is the item parameters plus the
# K-1 size logits, S11 is the whole of .step1_variance()'s matrix, and theta2
# is the outcome model on the packing the finish uses (.qn_pack_sm()).
#
# References
#   Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#     classes and external variables. Psychometrika, 83(4), 871-892.
#   Gong, G., & Samaniego, F. J. (1981). Pseudo maximum likelihood estimation:
#     Theory and applications. The Annals of Statistics, 9(4), 861-869.

# The blocks of the joint information matrix a two-step covariate fit needs,
# by central differences on the second-difference step `.fd_hessian_step`
# (the scores share its perturbed densities). Returns H22 (p2 x p2) and H12
# (p1 x p2), both as second derivatives of the log-likelihood (negate for
# information), and the n x p2 matrix of case-level scores in theta2. H12 is
# NULL when the measurement model has no unconstrained packing. `pack` and
# `unpack` give theta2; on a structural model that does not supply the class
# probabilities (a distal outcome) the measurement side carries the step-1
# class sizes and theta1 is .step1_pack()'s vector, sizes included.
.twostep_information <- function(model_state, X, Y, w, pack = .step1_pack_sm,
                                 unpack = .step1_unpack_sm) {
  mm <- model_state$mm
  sm <- model_state$sm
  th2 <- pack(sm)
  p2  <- length(th2)
  h2  <- .fd_hessian_step * pmax(1, abs(th2))

  sizes <- !.supplies_class_probs(sm)
  A_of  <- function(ms) {
    A <- log_likelihood(ms$mm, X)
    if (sizes) A <- sweep(A, 2, log(pmax(ms$weights, 1e-300)), "+")
    A
  }
  A0   <- A_of(model_state)
  B_at <- function(v) log_likelihood(unpack(sm, v), Y)
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

  th1 <- if (sizes) .step1_pack(model_state) else .step1_pack_mm(mm)
  if (is.null(th1) || !length(th1))
    return(list(H22 = H22, H12 = NULL, S2 = S2))
  un1 <- if (sizes) function(v) .step1_unpack(model_state, v) else
    function(v) list(mm = .step1_unpack_mm(mm, v))

  p1  <- length(th1)
  h1  <- .fd_hessian_step * pmax(1, abs(th1))
  H12 <- matrix(0, p1, p2)
  for (a in seq_len(p1)) {
    va <- th1; va[a] <- va[a] + h1[a]
    vb <- th1; vb[a] <- vb[a] - h1[a]
    Ap <- A_of(un1(va))
    Am <- A_of(un1(vb))
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
# @param w Case weights of step 2.
# @param se One of "corrected", "robust", "hessian".
# @param w1 Case weights of step 1: `w` except where cases missing a class
#   predictor are listwise deleted, which weigh zero in `w` only.
#
# @return A list with `V` (the K*D square covariance matrix, the anchored class
#   padded with zeros) and `method`, a label naming the estimator used.
.twostep_covariate_vcov <- function(model_state, X, Y, w, se = "corrected",
                                    w1 = w) {
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
  # theta2 maximises the likelihood plus the latent prior the class regression
  # is fitted with, so its information carries the prior's curvature too (the
  # note above .step3_scores(), R/step3_variance.R). th2 is .step1_pack_sm()'s
  # order, beta recentred on row K, which is the helper's.
  beta0 <- sweep(sm$parameters$beta, 2, sm$parameters$beta[K, ], "-")
  Hp   <- .covariate_prior_derivs(beta0, .covariate_design(sm, Y),
                                  .bayes_alpha(sm, "latent"), w)$hessian
  V2   <- pinv(-(info$H22 + Hp))
  V2   <- (V2 + t(V2)) / 2

  if (identical(se, "hessian"))
    return(list(V = pad(V2), method = "Two-step observed information (step 2 only)"))

  if (identical(se, "robust")) {
    survey <- isTRUE(sm$has_survey_design) && length(sm$strata) == nrow(X) &&
      length(sm$cluster) == nrow(X)
    # The scores carry w; a frequency weight's meat takes sqrt(w) instead
    # (.meat_weights()), and counts cases rather than rows.
    wtype <- model_state$weight_type %||% "sampling"
    S2 <- info$S2 * (.meat_weights(w, wtype) / pmax(w, 1e-300))
    D2 <- V2 %*% (.score_meat(S2, sm$strata, sm$cluster, survey,
                              n = .meat_n(w, wtype)) - Hp) %*% V2
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
                               "step-1 term unavailable for this measurement model"),
                step1_unavailable = TRUE))

  p1 <- length(th1)
  d1_method <- if (p1 <= .step1_hessian_max) "hessian" else "outer"
  if (d1_method == "outer")
    message(sprintf(paste("Step-1 variance for the two-step standard errors",
                          "uses the outer-product estimator: the measurement",
                          "model has %d parameters (limit %d for the numerical",
                          "Hessian)."), p1, .step1_hessian_max))
  S1  <- .step1_case_scores(model_state, X, th1, w1)
  d1  <- .step1_variance(model_state, X, th1, S1, w1, d1_method)
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

# The outcome models the distal two-step variance covers: a lone continuous
# outcome (class means and variances) and a lone categorical one (class
# logits). Outcomes with covariates keep the estimator their M-step attaches.
.twostep_distal_ok <- function(sm, Y) {
  if (inherits(sm, "distal_continuous")) return(TRUE)
  inherits(sm, "distal_categorical") && NCOL(Y) == 1L
}

# Variance-covariance matrix of a two-step distal outcome, on .qn_pack_sm()'s
# vector theta2, with the same three `se` choices as the covariate case. Step
# 2 maximises the joint likelihood plus the outcome's own prior
# (.qn_sm_prior(), zero for a continuous outcome), so its curvature joins H22.
# Under "corrected", G = V2 I12' is d theta2_hat / d theta1 (the implicit
# function of step 2's score equations), V = V2 + G S11 G', and `Vq` is the
# joint covariance of (size logits, theta2): the sizes are step 1's, so their
# block is S11's and their covariance with theta2 is G S11. Under "hessian" and
# "robust" step 1 is taken as known and the size block is zero.
.twostep_distal_vcov <- function(model_state, X, Y, w, se = "corrected", w1 = w) {
  K   <- model_state$n_components
  sm  <- model_state$sm
  pk  <- function(s) .qn_pack_sm(s)$par
  th2 <- pk(sm)
  p2  <- length(th2)
  info <- .twostep_information(model_state, X, Y, w, pk, .qn_unpack_sm)
  lp   <- .qn_sm_prior(sm, Y, w)
  Hp   <- .step1_fd_hessian(function(v) lp(.qn_unpack_sm(sm, v)), th2,
                            .fd_hessian_step)
  V2 <- pinv(-(info$H22 + Hp))
  V2 <- (V2 + t(V2)) / 2
  q  <- seq_len(K - 1L)
  joint <- function(V, Sqq = matrix(0, K - 1L, K - 1L),
                    C = matrix(0, p2, K - 1L))
    rbind(cbind(Sqq, t(C)), cbind(C, V))

  if (identical(se, "hessian"))
    return(list(V = V2, Vq = joint(V2),
                method = "Two-step observed information (step 2 only)"))

  if (identical(se, "robust")) {
    survey <- isTRUE(sm$has_survey_design) && length(sm$strata) == nrow(X) &&
      length(sm$cluster) == nrow(X)
    wtype <- model_state$weight_type %||% "sampling"
    S2 <- info$S2 * (.meat_weights(w, wtype) / pmax(w, 1e-300))
    D2 <- V2 %*% (.score_meat(S2, sm$strata, sm$cluster, survey,
                              n = .meat_n(w, wtype)) - Hp) %*% V2
    D2 <- (D2 + t(D2)) / 2
    return(list(V = D2, Vq = joint(D2),
                method = sprintf("Two-step sandwich (%s)",
                                 if (isTRUE(sm$has_survey_design))
                                   "survey-linearized" else "robust")))
  }

  th1 <- .step1_pack(model_state)
  if (is.null(info$H12) || is.null(th1) || !length(th1))
    return(list(V = V2, Vq = joint(V2),
                method = paste("Two-step observed information (step 2 only);",
                               "step-1 term unavailable for this measurement model"),
                step1_unavailable = TRUE))

  p1 <- length(th1)
  d1_method <- if (p1 <= .step1_hessian_max) "hessian" else "outer"
  if (d1_method == "outer")
    message(sprintf(paste("Step-1 variance for the two-step standard errors",
                          "uses the outer-product estimator: the measurement",
                          "model has %d parameters (limit %d for the numerical",
                          "Hessian)."), p1, .step1_hessian_max))
  S1  <- .step1_case_scores(model_state, X, th1, w1)
  d1  <- .step1_variance(model_state, X, th1, S1, w1, d1_method)
  S11 <- d1$V
  G   <- V2 %*% t(info$H12)
  V   <- V2 + G %*% S11 %*% t(G)
  V   <- (V + t(V)) / 2
  iq  <- p1 - (K - 1L) + q
  list(V = V,
       Vq = joint(V, S11[iq, iq, drop = FALSE], G %*% S11[, iq, drop = FALSE]),
       method = sprintf("Two-step pseudo-ML (Bakk and Kuha, 2018; %s step 1)%s",
                        d1$method,
                        if (isTRUE(d1$fallback))
                          "; step-1 Hessian was not positive definite" else ""))
}

# Write a two-step distal variance where the step-3 readers look for one
# (.attach_step3_distal_vcov(), R/corrections.R): a continuous outcome's
# `Sigma_mu`, `ses` and the (sizes, means, variances) matrix `V_pmv` that
# outcome_contrasts() takes its effect sizes from; a categorical outcome's
# `hessian`, on the same as.vector(beta_pooled) layout.
.store_twostep_distal <- function(sm, res, weights) {
  K <- sm$n_components
  V <- res$V
  if (inherits(sm, "distal_continuous")) {
    m <- seq_len(K)
    sm$parameters$Sigma_mu <- V[m, m, drop = FALSE]
    sm$parameters$ses      <- matrix(sqrt(pmax(diag(V)[m], 0)), ncol = 1)
    pk  <- pmax(weights, 1e-12)
    par <- c(log(pk[-K] / pk[K]), .qn_pack_sm(sm)$par)
    ia  <- seq_len(K - 1L)
    to_pmv <- function(th) {
      a <- c(th[ia], 0); p <- exp(a - max(a)); p <- p / sum(p)
      u <- .qn_unpack_sm(sm, th[-ia])$parameters
      c(p, u$means[, 1L], rep_len(u$covariances[, 1L], K))
    }
    h <- 1e-4 * pmax(1, abs(par))
    J <- vapply(seq_along(par), function(j) {
      a <- par; a[j] <- a[j] + h[j]
      b <- par; b[j] <- b[j] - h[j]
      (to_pmv(a) - to_pmv(b)) / (2 * h[j])
    }, numeric(3L * K))
    sm$parameters$V_pmv    <- J %*% res$Vq %*% t(J)
    sm$parameters$pk_step3 <- to_pmv(par)[seq_len(K)]
  } else {
    sm$parameters$hessian <- -pinv(V)
  }
  sm$parameters$V_method <- res$method
  sm
}

# Compute and attach the two-step variance to the structural sub-model, the
# way .attach_step3_covariate_vcov() does for the third step: a class
# regression's on `V_robust`, a lone distal outcome's through
# .store_twostep_distal(). Any other structural model keeps the variance its
# M-step attached.
.attach_twostep_covariate_vcov <- function(model_state, X, Y, w, se = "corrected",
                                           w1 = w) {
  if (identical(se, "none")) return(model_state)
  if (.twostep_distal_ok(model_state$sm, Y)) {
    res <- tryCatch(
      .twostep_distal_vcov(model_state, X, Y, w, se = se, w1 = w1),
      error = function(e) {
        warning(sprintf(
          "Two-step outcome standard errors fell back to the M-step's: %s",
          conditionMessage(e)), call. = FALSE)
        NULL
      })
    if (is.null(res)) return(model_state)
    if (identical(se, "corrected") && isTRUE(res$step1_unavailable))
      .warn_step1_unavailable(model_state$mm,
                              reported = "step 2's observed information",
                              instead  = "hessian")
    model_state$sm <- .store_twostep_distal(model_state$sm, res,
                                            model_state$weights)
    return(model_state)
  }
  if (!inherits(model_state$sm, "covariate")) return(model_state)

  res <- tryCatch(
    .twostep_covariate_vcov(model_state, X, Y, w, se = se, w1 = w1),
    error = function(e) {
      warning(sprintf(
        "Two-step covariate standard errors fell back to the Q-function Hessian: %s",
        conditionMessage(e)), call. = FALSE)
      NULL
    })
  if (is.null(res)) return(model_state)
  if (identical(se, "corrected") && isTRUE(res$step1_unavailable))
    .warn_step1_unavailable(model_state$mm,
                            reported = "step 2's observed information",
                            instead  = "hessian")

  model_state$sm$parameters$V_robust <- res$V
  model_state$sm$parameters$V_method <- res$method
  model_state
}

# ==============================================================================
# The same variance for a latent transition model
# ==============================================================================
#
# fit_lta(n_steps = 2) (R/lta.R) freezes the measurement block at its step-one
# estimate and maximises the initial-status and transition regressions on the
# full likelihood, so it is the estimator above in a different model and it
# takes the same variance. What changes is the packing: an lta_model is not
# covered by .step1_pack(), and does not need to be -- .lta_par_layout() /
# .lta_par_pack() / .lta_ll_case() (R/lta_core.R) already describe the whole
# free parameter vector and give a case-level log-likelihood on it, and every
# block of that layout carries a `kind`, so the (theta1, theta2) partition is
# .lta_par_split()'s filter on that tag rather than a second packing.
#
# One difference from the cross-sectional case is worth stating. There the
# joint log-likelihood factors into a measurement piece and a structural piece
# that can be perturbed separately and recombined, which is why
# .twostep_information() re-evaluates only one of them at a time. The forward-
# backward recursion does not factor that way: a latent transition model's
# case-level log-likelihood mixes the two blocks at every occasion, so each
# perturbed value costs a full pass. The blocks are therefore differenced
# directly, and only the two that equation 5 uses: I22 (p2 x p2) and I12
# (p1 x p2). I11 is never formed. On a real model p1 is several times p2 and
# that corner is most of the cost of a full Hessian, which the robust branch of
# .lta_standard_errors() pays because a sandwich needs it and this does not.

# I22 and I12 of the joint log-likelihood at the fitted point, by the same
# four-point central difference .step1_fd_hessian() uses, on its relative
# step `.fd_hessian_step` (R/step3_variance.R, which says why it is 2e-4). The
# prior's curvature, which every latent transition variance adds, is
# differenced on the same step: the sandwich magnifies any noise in it.
# Returned as second derivatives of the log-likelihood; negate for
# information. `ll` is the objective on the whole vector; the default is the
# joint log-likelihood of `state`, and the three-step passes its own, in which
# the step-1 coordinates reach step 3 through the classification error.
.lta_twostep_information <- function(state, X, w, layout, par, idx1, idx2,
                                     ll = function(v)
                                       sum(w * .lta_ll_case(state, X, v, layout))) {
  h  <- .fd_hessian_step * pmax(1, abs(par))
  bump <- function(i, s) { v <- par; v[i] <- v[i] + s * h[i]; v }
  two  <- function(i, si, j, sj) {
    v <- par; v[i] <- v[i] + si * h[i]; v[j] <- v[j] + sj * h[j]; ll(v)
  }
  cross <- function(i, j)
    (two(i, 1, j, 1) - two(i, 1, j, -1) - two(i, -1, j, 1) + two(i, -1, j, -1)) /
      (4 * h[i] * h[j])

  p2  <- length(idx2)
  H22 <- matrix(0, p2, p2)
  for (a in seq_len(p2)) for (b in a:p2) {
    if (a == b) {
      # The diagonal is the same four-point formula with i = j, which collapses
      # to the ordinary second difference; writing it out avoids evaluating the
      # unperturbed point twice as a mixed term.
      H22[a, a] <- (ll(bump(idx2[a], 2)) - 2 * ll(par) + ll(bump(idx2[a], -2))) /
        (4 * h[idx2[a]]^2)
    } else {
      H22[a, b] <- H22[b, a] <- cross(idx2[a], idx2[b])
    }
  }

  p1  <- length(idx1)
  H12 <- matrix(0, p1, p2)
  for (a in seq_len(p1)) for (b in seq_len(p2))
    H12[a, b] <- cross(idx1[a], idx2[b])

  list(H22 = H22, H12 = H12)
}

# The prior's Hessian at `par`, by central differences of its closed-form
# gradient `grad` (.lta_penalty()'s, which runs no forward-backward pass): 2p
# evaluations where differencing the prior `lp` itself costs 2p(p + 1), 21 s
# of the ENPOL joint model's 95 at p = 147. The gradient is used only where it
# agrees with `lp`'s own differences at `par`; elsewhere, or where it cannot
# be written down (NULL), `lp` is differenced as before.
.lta_prior_hessian <- function(lp, grad, par) {
  g0 <- tryCatch(grad(par), error = function(e) NULL)
  ok <- length(g0) == length(par) && all(is.finite(g0))
  if (ok) {
    fd <- .qn_fd_grad2(lp, par)$g
    ok <- max(abs(g0 - fd)) <= 1e-5 * (1 + max(abs(fd)))
  }
  if (!ok) return(.step1_fd_hessian(lp, par, .fd_hessian_step))
  H <- vapply(seq_along(par), function(j) {
    h <- 1e-5 * max(1, abs(par[j]))
    e <- numeric(length(par)); e[j] <- h
    (grad(par + e) - grad(par - e)) / (2 * h)
  }, numeric(length(par)))
  (H + t(H)) / 2
}

# Columns `cols` of the log-likelihood's Hessian on a packed vector `par`, by
# central differences of its analytic gradient (the case scores summed under
# the case weights): 2 x length(cols) score passes, where differencing the
# log-likelihood itself costs four passes per entry. This is
# .lta_penalised_hessian()'s construction (R/lta.R), on its step, restricted
# to the columns a caller needs. `par` is the free vector; the packed vector
# of `layout` is `par[src]`, so `src` also folds the score columns back onto
# `par` where several of them share one free coordinate (the shared
# prevalence vector of .qn_lta_reduction()), as the finish folds them.
#
# The passes run on one row per distinct (responses, covariates) pattern
# (.qn_lta_pattern_rows()), unpacking straight onto it, kept only where that
# table reproduces the full sample's log-likelihood at `par`; the scores are
# sums over cases that merely repeat, so the derivatives are the full
# sample's to the rounding of the sums. On the ENPOL joint models that is
# 15,347 rows in 60,784. The columns are independent, so they are cut into
# one batch per worker on `n_cores`; each batch carries the pattern table
# only, and the answer is the same at every `n_cores`.
#
# `occasion`, for the whole Hessian only (`cols` every coordinate), names the
# one occasion each coordinate's data comes from, 0 for a coordinate that
# reaches more than one. Coordinates of different occasions of a fit whose
# occasions are independent have no cross curvature, so one pass perturbs one
# coordinate of each occasion at once and each occasion's rows of the
# difference are that coordinate's column; a 0 coordinate gets a pass of its
# own, and its rows in the others' columns are read off its column by
# symmetry. The ENPOL joint models' step 1 is 83 + 64 coordinates: 83 passes
# in place of 147 (.lta_shared_occasions()).
#
# NULL where the scores do not describe every coordinate (the caller then
# differences the log-likelihood), and where a frozen fit's emission tables
# are cached on the state, since a perturbed measurement parameter would not
# move them.
.lta_score_hessian_cols <- function(state, X, layout, par, cols,
                                    src = seq_along(par), shared_ok = FALSE,
                                    n_cores = 1L, occasion = NULL) {
  if (!is.null(state$logB_fixed)) return(NULL)
  base <- state; w <- state$weights_vec
  pr <- .qn_lta_pattern_rows(state, X)
  if (!is.null(pr)) {
    rows_st <- pr$shrink(state)
    full <- sum(w * .lta_ll_case_state(
      .lta_par_unpack(par[src], state, layout, average = FALSE), X))
    rows <- sum(pr$w * .lta_ll_case_state(
      .lta_par_unpack(par[src], rows_st, layout, average = FALSE), pr$X))
    if (isTRUE(abs(full - rows) <= 1e-8 * (1 + abs(full)))) {
      base <- rows_st; X <- pr$X; w <- pr$w
    }
  }
  # What each batch is sent. A fitted state carries its call, with the data
  # inline when it came through do.call(), and the data again; neither is read
  # by the scores, and on the ENPOL joint models they were 126 of the 129 MB
  # a batch would otherwise carry to every worker.
  base$call <- NULL; base$data <- NULL
  if (is.null(.lta_score_grad(par, base, X, w, layout, src, shared_ok)))
    return(NULL)
  if (!identical(as.integer(cols), seq_along(par)) ||
      length(occasion) != length(par)) occasion <- NULL
  probes <- if (is.null(occasion)) as.list(cols) else {
    by <- split(cols[occasion > 0L], occasion[occasion > 0L])
    c(as.list(cols[occasion == 0L]),
      lapply(seq_len(max(0L, lengths(by))), function(i)
        unlist(lapply(by, function(v) if (i <= length(v)) v[i]),
               use.names = FALSE)))
  }
  # The pool is sized by `n_cores`, as the fit's own search sized it, so the
  # cached workers are reused; only the number of batches follows the passes.
  n_b    <- max(1L, min(suppressWarnings(as.integer(n_cores)), length(probes)),
                na.rm = TRUE)
  batch  <- rep(seq_len(n_b), length.out = length(probes))
  pieces <- split(probes, batch)
  out <- .par_lapply(pieces, .lta_score_hessian_piece, n_cores = n_cores,
                     base = base, Y = X, w = w, layout = layout, par = par,
                     src = src, shared_ok = shared_ok)
  if (any(vapply(out, is.null, logical(1)))) return(NULL)
  D <- do.call(cbind, out)[, order(order(batch)), drop = FALSE]

  h <- 1e-5 * pmax(1, abs(par))
  M <- matrix(0, length(par), length(cols))
  for (s in seq_along(probes)) for (j in probes[[s]]) {
    rows <- if (is.null(occasion) || occasion[j] == 0L) seq_along(par) else
      which(occasion == occasion[j])
    M[rows, match(j, cols)] <- D[rows, s] / h[j]
  }
  z <- which(occasion == 0L)
  if (length(z)) {
    k <- which(occasion > 0L)
    M[z, k] <- t(M[k, z, drop = FALSE])
  }
  M
}

# The summed scores at `par` (packed as `par[src]`, folded back by `src`), or
# NULL where the score blocks do not describe every coordinate.
.lta_score_grad <- function(par, base, X, w, layout, src, shared_ok) {
  st <- .lta_par_unpack(par[src], base, layout, average = FALSE)
  sc <- .lta_score_matrix(st, X, shared_ok = shared_ok, weights = w)
  if (is.null(sc) || ncol(sc$S) != length(src)) return(NULL)
  as.vector(rowsum(sc$S[1L, ], src))
}

# One worker's batch of .lta_score_hessian_cols(): for each probe (a set of
# coordinates perturbed together, each on .lta_penalised_hessian()'s
# relative step), half the central difference of the summed scores. The data
# is `Y`, not `X`, which .par_lapply() would take for the list it maps over.
.lta_score_hessian_piece <- function(probes, base, Y, w, layout, par, src,
                                     shared_ok) {
  D <- matrix(0, length(par), length(probes))
  for (i in seq_along(probes)) {
    j <- probes[[i]]
    e <- numeric(length(par)); e[j] <- 1e-5 * pmax(1, abs(par[j]))
    gp <- .lta_score_grad(par + e, base, Y, w, layout, src, shared_ok)
    gm <- .lta_score_grad(par - e, base, Y, w, layout, src, shared_ok)
    if (is.null(gp) || is.null(gm)) return(NULL)
    D[, i] <- (gp - gm) / 2
  }
  D
}

# The occasion each coordinate of a transition-free fit's reduced vector
# reads its data from, for .lta_score_hessian_cols()'s paired passes: the
# initial-status blocks the first, the transition blocks the second, an item
# block the one occasion the item is observed at (the joint model's items are
# each one variable's, the other occasion padded with NA), and 0 for anything
# else. NULL unless the occasions are independent: two of them, one class,
# no random intercept. Item slopes on ordinal items run on the random
# intercept's code with one node at zero and every loading fixed at zero
# (R/lta_dif.R), which links nothing, so that one counts as none.
.lta_shared_occasions <- function(red, X) {
  st <- red$state
  ri <- st$ri
  linked <- !is.null(ri) && (NROW(ri$Dnode) != 1L || any(ri$Dnode != 0) ||
                               any(ri$L != 0) || isTRUE(ri$loading_free))
  if (!identical(as.integer(st$n_times), 2L) || (st$n_classes %||% 1L) != 1L ||
      linked) return(NULL)
  J  <- st$n_items
  at <- lapply(seq_len(J), function(j)
    which(vapply(1:2, function(t) any(!is.na(X[, .time_block_cols(t, J)[j]])),
                 logical(1))))
  g <- unlist(lapply(red$layout, function(b) {
    o <- if (b$kind %in% c("delta", "delta_beta")) 1L
    else if (b$kind %in% c("tau", "tau_beta")) 2L
    else if (!is.null(b$j) && length(b$j) == 1L && length(at[[b$j]]) == 1L)
      at[[b$j]]
    else 0L
    rep(o, b$len)
  }))
  if (length(g) != length(red$src)) return(NULL)
  # A reduced coordinate that several packed columns share keeps an occasion
  # only if they all agree.
  vapply(split(g, red$src), function(v) if (all(v == v[1])) v[1] else 0L,
         integer(1))
}

# The curvature of a transition-free fit's objective on the vector it actually
# has. That model ties every origin row of an occasion to one prevalence
# vector, which the ordinary packing does not describe, so it is read through
# .qn_lta_reduction(): the log-likelihood and the prior are differenced on the
# reduced vector `red$r0`, and `V` is the inverse of that penalised curvature.
# Taken over the response patterns unless the fit reads a covariate (class
# regressions or item slopes), which makes every case its own pattern. `Hp` is
# the prior's share, kept for a sandwich's meat. Used for the step-1 variance
# of the joint model's two-step and for the standard errors of the joint
# model with no association.
#
# The curvature depends on step 1 alone, and the step-2 models fitted on one
# step 1 (the ENPOL joint models m0 ... m2x share theirs) would each pay for
# it again: 166 of m0's 194 score passes. The last one is kept and handed back
# when the same reduced model, data, weights and smoothing come in again,
# compared whole by identical(). Left out of the comparison are only what the
# curvature never reads -- the call, the reported fit summaries (`metrics`,
# `n_params`, `n_iter`) and the joint model's labels, which a step 1 handed
# over to a later fit recounts or carries and a fresh one does not --
# and the data and per-case results the reduced model's parameters and `X`
# already determine. A field missed here costs a recomputation, never a wrong
# answer. Data above 1e7 cells is not kept, as .one_hot_memo() does not keep
# it.
.lta_curvature_cache <- new.env(parent = emptyenv())

.lta_shared_curvature <- function(state, X, alpha, red = .qn_lta_reduction(state),
                                  n_cores = 1L) {
  if (is.null(red)) return(NULL)
  s0  <- unclass(red$state)
  skip <- c("call", "metrics", "n_params", "n_iter", "n_x", "n_y",
            "n_items_x", "n_items_y", "joint_stacked", "joint_key",
            "data", "gamma", "ll_case")
  key <- list(state = s0[sort(setdiff(names(s0), skip))],
              layout = red$layout, r0 = red$r0, src = red$src, X = X,
              w = state$weights_vec, alpha = alpha)
  hit <- .lta_curvature_cache$last
  if (!is.null(hit) && identical(hit$key, key)) {
    out <- hit$out; out$red <- red
    return(out)
  }
  out <- .lta_shared_curvature_fit(state, X, alpha, red, n_cores)
  if (length(X) <= 1e7) .lta_curvature_cache$last <- list(key = key, out = out)
  out
}

.lta_shared_curvature_fit <- function(state, X, alpha, red, n_cores) {
  layout <- red$layout
  s1 <- red$state
  per_case <- !is.null(s1$Z_delta) || !is.null(s1$Z_tau) || !is.null(s1$dif)
  if (per_case) {
    Xc <- X; wc <- state$weights_vec
  } else {
    pat <- .pattern_index(X)
    Xc <- X[pat$rep_row, , drop = FALSE]
    wc <- as.vector(rowsum(state$weights_vec, pat$idx))
  }
  s1$weights_vec <- wc
  ll <- function(r) sum(wc * .lta_ll_case(s1, Xc, r[red$src], layout))
  lp <- function(r)
    .lta_log_prior(.lta_par_unpack(r[red$src], s1, layout, average = FALSE),
                   Xc, alpha)
  lg <- function(r)
    as.vector(rowsum(.lta_penalty(.lta_par_unpack(r[red$src], s1, layout,
                                                  average = FALSE),
                                  Xc, layout, alpha)$gradient, red$src))
  Hp <- if (.lta_prior_is_zero(s1, alpha) || is.na(lp(red$r0))) 0 else
    .lta_prior_hessian(lp, lg, red$r0)
  # The log-likelihood's curvature from the differenced scores, folded onto
  # the reduced vector: 2 x length(r0) score passes in place of the
  # log-likelihood differenced over every pair. The likelihood route stays
  # for a fit the scores do not cover.
  Hd <- tryCatch(.lta_score_hessian_cols(red$state, X, layout, red$r0,
                                         seq_along(red$r0), red$src,
                                         shared_ok = TRUE, n_cores = n_cores,
                                         occasion = .lta_shared_occasions(red, X)),
                 error = function(e) NULL)
  Hd <- if (!is.null(Hd) && all(is.finite(Hd))) (Hd + t(Hd)) / 2 else
    .step1_fd_hessian(ll, red$r0)
  V <- .psd_pinv(-(Hd + Hp))
  list(V = (V + t(V)) / 2, Hp = Hp, red = red)
}

# Standard errors of a transition-free fit with no transition regression: the
# joint model with no association (fit_joint_classes(association = "none")).
# The variance on the reduced vector is expanded back onto the ordinary
# packing, `V[src, src]`, so every origin row of the shared prevalence vector
# carries the same variance and every block the score matrix names reads the
# columns it always reads. `robust` is the sandwich on the case scores folded
# onto the reduced vector, the prior's curvature taken out of the meat as in
# .lta_standard_errors(). With a transition regression the origin-free design
# packs as it is and .lta_standard_errors() serves the fit directly.
.lta_shared_standard_errors <- function(state, X, robust = FALSE, alpha = 0,
                                        n_cores = 1L) {
  cv <- tryCatch(.lta_shared_curvature(state, X, alpha, n_cores = n_cores),
                 error = function(e) NULL)
  if (is.null(cv) || !all(is.finite(cv$V))) return(NULL)
  red <- cv$red
  sc  <- .lta_score_matrix(red$state, X, shared_ok = TRUE)
  if (is.null(sc) || ncol(sc$S) != length(red$src)) return(NULL)
  Vr <- cv$V
  if (isTRUE(robust)) {
    w  <- state$weights_vec
    Sr <- t(rowsum(t(sc$S), red$src))
    n  <- state$n_eff %||% length(w)
    B  <- (n / max(n - 1, 1)) * (t(Sr) %*% sweep(Sr, 1, w, "*"))
    Vr <- Vr %*% (B - cv$Hp) %*% Vr
    Vr <- (Vr + t(Vr)) / 2
  }
  V <- Vr[red$src, red$src, drop = FALSE]
  list(vcov = V, blocks = sc$blocks, prob_se = .lta_prob_se(sc$blocks, V),
       loading_se = NULL, conditional = sc$conditional, design_based = FALSE,
       robust = isTRUE(robust), twostep = FALSE, twostep_V2 = NULL,
       method = if (isTRUE(robust))
         "Robust sandwich around the observed information"
       else "Observed information")
}

# Pseudo-maximum-likelihood variance of a two-step latent transition fit.
#
# @param state The step-two fit: measurement block at its step-one values,
#   structural block at the joint maximum.
# @param X Indicator matrix both steps were fitted to.
# @param step1 The step-one fit, as carried on the object as `$step1`.
# @param alpha The `smoothing` both steps were fitted with.
#
# @return NULL when the model is outside the packing, or a list of `V` (the
#   structural block's covariance matrix), `V2` (its step-two-only part), `S11`
#   (the step-one variance actually used), the two index vectors and a `method`
#   label. Returning NULL rather than signalling is deliberate: a model the
#   packing cannot describe keeps the estimator it already had, the way
#   .lta_standard_errors()'s robust branch declines.
.lta_twostep_vcov <- function(state, X, step1, alpha = 0, n_cores = 1L) {
  if (is.null(step1) || !inherits(step1, "lta_model")) return(NULL)
  # A transition-free step 1 (the joint model's, fit_joint_classes()) ties
  # every origin row to one prevalence vector, which the ordinary packing does
  # not describe; it is read through the reduction .lta_threestep_vcov() uses.
  indep <- isTRUE(step1$tau_independent)
  if (!.lta_par_packable(state) || (!indep && !.lta_par_packable(step1)))
    return(NULL)

  layout <- .lta_par_layout(state)
  par    <- .lta_par_pack(state, layout)
  split  <- .lta_par_split(layout)
  idx1   <- split$measurement
  idx2   <- split$structural
  if (!length(idx1) || !length(idx2)) return(NULL)

  red <- if (indep) .qn_lta_reduction(step1)
  if (indep && is.null(red)) return(NULL)
  layout1 <- if (indep) red$layout else .lta_par_layout(step1)
  par1    <- if (indep) red$r0[red$src] else .lta_par_pack(step1, layout1)
  idx1_s1 <- .lta_par_split(layout1)$measurement

  # The freeze invariant, used here as an alignment proof. Both vectors come
  # from the same layout code on the same model shape, so their measurement
  # halves are the same blocks in the same order; if step two really did leave
  # the measurement block alone, the two halves are also numerically identical.
  # Anything else -- a status relabelling between the two fits, a layout that
  # does not correspond -- shows up here, and the honest response is to decline
  # rather than to propagate a variance read off the wrong coordinates.
  if (length(idx1_s1) != length(idx1)) return(NULL)
  if (max(abs(par[idx1] - par1[idx1_s1])) > 1e-6) return(NULL)

  # Sigma_11, the sampling variance of the step-one estimate, restricted to the
  # measurement parameters: the corresponding block of the INVERSE step-one
  # information, not the inverse of the block. The step-one fit carries it
  # already when it was asked for standard errors; otherwise it is computed
  # here, once.
  # The transition-free step 1 reports no standard errors, so its variance is
  # the inverse penalised curvature on the reduced vector, with the measurement
  # coordinates read back one-to-one. It is taken over the response patterns
  # unless step 1 reads a covariate (item slopes, and the class regressions
  # that come with them), which makes every case its own pattern.
  if (indep) {
    S1r <- .lta_shared_curvature(step1, as.matrix(step1$data), alpha, red,
                                 n_cores)$V
    j   <- red$src[idx1_s1]
    S11 <- S1r[j, j, drop = FALSE]
  } else {
    se1 <- step1$se
    if (is.null(se1) || is.null(se1$vcov))
      se1 <- tryCatch(.lta_standard_errors(step1, X, alpha = alpha),
                      error = function(e) NULL)
    if (is.null(se1) || is.null(se1$vcov)) return(NULL)
    if (nrow(se1$vcov) != length(par1)) return(NULL)
    S11 <- se1$vcov[idx1_s1, idx1_s1, drop = FALSE]
  }

  # The columns along the structural coordinates of the Hessian, from the
  # differenced scores: 2 x p2 score passes, against 4 x p1 x p2 passes of
  # the log-likelihood for H12 alone, which ran past 80 minutes on the ENPOL
  # joint models (n 60,784, p1 133). The likelihood route stays for a model
  # the scores do not cover.
  M <- tryCatch(.lta_score_hessian_cols(state, X, layout, par, idx2,
                                        n_cores = n_cores),
                error = function(e) NULL)
  info <- if (!is.null(M) && all(is.finite(M))) {
    H22 <- M[idx2, , drop = FALSE]
    list(H22 = (H22 + t(H22)) / 2, H12 = M[idx1, , drop = FALSE])
  } else .lta_twostep_information(state, X, state$weights_vec, layout, par,
                                  idx1, idx2)
  # Step 2 maximises the log-likelihood plus the structural prior, so its
  # information is that objective's: the prior's curvature joins H22. The
  # prior reads the structural vector and the fixed items only, so the cross
  # derivative H12 is the likelihood's alone. S11 above is step 1's own
  # variance, the penalised curvature of the objective step 1 maximised.
  marg <- if (is.null(state$ri))
    tryCatch(.lta_prior_marginals(state, X), error = function(e) NULL)
  lp2 <- function(v) {
    p <- par; p[idx2] <- v
    .lta_log_prior(.lta_par_unpack(p, state, layout, average = FALSE), X,
                   alpha, marg)
  }
  lg2 <- function(v) {
    p <- par; p[idx2] <- v
    .lta_penalty(.lta_par_unpack(p, state, layout, average = FALSE), X, layout,
                 alpha)$gradient[idx2]
  }
  Hp22 <- if (.lta_prior_is_zero(state, alpha) || is.na(lp2(par[idx2]))) 0 else
    .lta_prior_hessian(lp2, lg2, par[idx2])
  V2 <- .psd_pinv(-(info$H22 + Hp22))
  V2 <- (V2 + t(V2)) / 2
  V1 <- V2 %*% t(info$H12) %*% S11 %*% info$H12 %*% V2
  V  <- V2 + V1
  V  <- (V + t(V)) / 2
  if (!all(is.finite(V))) return(NULL)

  list(V = V, V2 = V2, S11 = S11, idx1 = idx1, idx2 = idx2,
       method = "Two-step pseudo-ML (Bakk and Kuha, 2018)")
}
