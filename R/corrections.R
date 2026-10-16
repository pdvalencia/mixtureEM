# ==============================================================================
# S3 3-Step Corrections (BCH and ML)
# ==============================================================================

get_modal_resp <- function(resp) {
  n <- nrow(resp)
  K <- ncol(resp)
  modal <- matrix(0, nrow = n, ncol = K)
  assigned_classes <- max.col(resp, ties.method = "first")
  # One-hot the modal class by index rather than a row loop. Also correct for
  # an empty `resp`, where `1:n` would have run over rows 1 and 0.
  modal[cbind(seq_len(n), assigned_classes)] <- 1
  return(modal)
}

# The step-two classification-error matrix of a three-step estimator.
#
# Step one fits a measurement model and step two turns its posteriors into an
# assigned-class variable W. That assignment is imperfect, and the third step
# needs to know how imperfect: `D[r, c] = P(W = r | true class = c)`, estimated
# as the posterior mass of class `c` that sits in cases assigned to `r`,
#
#     D[r, c] = sum_i A[i, r] resp[i, c] w_i / sum_i resp[i, c] w_i ,
#
# with `A` the assignment indicator -- the modal one-hot under
# `assignment = "modal"`, the posterior itself under "proportional". As in
# fit_bch() and fit_ml(), the *true* class is never observed, so its column
# weights are always the posteriors; only `A` changes with the rule. Rows are
# the assigned class and COLUMNS ARE THE TRUE CLASS, so columns sum to one.
#
#   ORIENTATION. This is the transpose of fit_ml()'s `C_row_norm` below, whose
#   rows are the true class. The two matrices are built from the same numbers
#   and neither is symmetric, so handing one where the other is wanted is
#   silent and wrong. Transpose deliberately at any point where they meet.
#
# The third step wants `D` as fixed multinomial logits with the last assigned
# category as the reference, `log(D[r, c] / D[K, c])`. A cell that is exactly
# zero makes that logit undefined, so `zero_floor` regularises it. Any value
# is arbitrary: the default 1e-6 is small enough to leave every other cell
# unchanged at the precision classification probabilities are reported to,
# while Nylund-Gibson, Grimm, Quirk and Furlong (2014) state 1e-4 in print.
# The argument is exposed rather than buried so that the choice, and the fact
# that it is a choice, is visible in the fitted object.
#
# `assigned`, when given, is the modal assignment to use instead of the one
# `resp` implies: an integer vector of assigned classes, one per row. The
# three-step variance differentiates `D` with respect to step 1's estimates,
# and the derivative of a modal table holds the assignments where they are --
# a small change to the posteriors must not move a case across a boundary, or
# the difference quotient measures the jump rather than the slope.
#
# Returns `D` (K x K) and `logits` (K x K, last row zero by construction).
#
# Reference
#   Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
#     transition mixture model using the three-step specification. Structural
#     Equation Modeling, 21(3), 439-454.
.classification_error <- function(resp,
                                  assignment = c("proportional", "modal"),
                                  weights    = NULL,
                                  zero_floor = 1e-6,
                                  assigned   = NULL) {

  assignment <- match.arg(assignment)
  resp <- as.matrix(resp)
  K    <- ncol(resp)

  if (nrow(resp) == 0L) stop("`resp` has no rows.")
  if (!is.finite(zero_floor) || zero_floor <= 0 || zero_floor >= 1 / K) {
    stop("`zero_floor` must be a small positive probability.")
  }

  if (is.null(weights)) weights <- rep(1, nrow(resp))
  if (length(weights) != nrow(resp)) {
    stop("`weights` must have one entry per row of `resp`.")
  }

  A  <- if (assignment == "proportional") resp else if (is.null(assigned))
    get_modal_resp(resp) else diag(K)[assigned, , drop = FALSE]
  RW <- resp * weights

  # Rows the assigned class, columns the true class. See ORIENTATION above.
  D  <- t(A) %*% RW
  Nk <- colSums(RW)

  # A class carrying no posterior mass at all leaves its column 0/0. It cannot
  # be estimated from these data, so say nothing about it rather than emit NaN.
  empty <- Nk <= 0
  D[, !empty] <- sweep(D[, !empty, drop = FALSE], 2, Nk[!empty], "/")
  D[, empty]  <- 1 / K

  # Floor, then restore the column sums the floor disturbed. Every column sums
  # to one before this line, so a column with no floored cell is unchanged.
  D <- pmax(D, zero_floor)
  D <- sweep(D, 2, colSums(D), "/")

  dimnames(D) <- list(assigned = seq_len(K), class = seq_len(K))

  logits <- log(sweep(D, 2, D[K, ], "/"))
  dimnames(logits) <- dimnames(D)

  list(D = D, logits = logits)
}

# Apply the BCH Correction
#
# Implementation follows the proportional BCH correction of Vermunt (2010)
# and Bakk, Tekle & Vermunt (2013):
#
#   Step 1: Build the classification error matrix C from proportional assignment.
#           C[k,j] = sum_i resp[i,k] * resp[i,j] / sum_i resp[i,j]
#           Columns index the true class and sum to 1.
#
#   Step 2: BCH weight matrix D = t(inv(C))
#
#   Step 3: Apply to proportional (posterior) weights:
#           bch_resp = resp %*% D
#
#   Negative weights are retained; they are a mathematically necessary
#   feature of the BCH correction for poorly-separated classes (Bakk et al.,
#   2013) and must not be clipped to zero.
#
# The BCH weights for assigned-class variable `A` and step-1 posteriors `resp`.
# The classification error matrix C[j, k] = P(assigned approx j | true = k)
# has columns indexing the true class, each summing to 1. The true class is
# never observed, so its column weights are always the posteriors; only the
# assigned variable A changes with the rule. Writing it as A crossed with
# `resp` rather than A crossed with itself matters: a modal A has a single 1
# per row, so t(A) %*% A is diagonal and would make C the identity, turning
# the correction into a no-op. Under proportional assignment A *is* `resp` and
# the two expressions are the same matrix. The weights are A %*% t(C^{-1});
# negative weights are retained intentionally. Kept in one place because the
# step-1 term of the joint variance (R/step3_joint.R) differentiates them.
.bch_weights <- function(resp, A, weights) {
  C <- t(A) %*% (resp * weights)
  C <- sweep(C, 2, colSums(resp * weights), "/")
  A %*% t(pinv(C))
}

fit_bch <- function(model_state, X, Y, assignment = "proportional",
                    max_iter = 1000, se = "robust") {

  if (inherits(model_state$sm, c("covariate", "distal_regression", "distal_pooled"))) {
    warning(paste(
      "BCH correction is not recommended for covariates or categorical distal outcomes.",
      "Vermunt (2010) and Bakk et al. (2013) show ML correction is preferred in this case,",
      "Consider correction = 'ML' instead."
    ))
  }

  # Complete missing covariates under the class-invariant Gaussian marginal
  # (endogenous-constrained-x; Sterba, 2014) so all cases are retained; a
  # missing distal outcome stays NA and is masked by the structural likelihood.
  Y <- complete_structural_covariates(model_state$sm, Y)

  weights <- model_state$sample_weights
  e_res   <- e_step(model_state, X, NULL)
  resp    <- exp(e_res$log_resp)          # n x K posterior probabilities

  # The assignment rule decides what the classification-error matrix is a table
  # of. Under proportional assignment a case contributes its posterior
  # probability to every class; under modal it contributes 1 to its most likely
  # class and 0 elsewhere (Bolck, Croon & Hagenaars, 2004; Vermunt, 2010). Both
  # are consistent; they differ in finite samples because the modal table sits
  # closer to the identity, so inverting it applies a smaller correction.
  A <- if (assignment == "modal") get_modal_resp(resp) else resp

  bch_resp <- .bch_weights(resp, A, weights)

  # A nested structural model (class predictors with distal outcomes) is not
  # always solved by one M-step -- a block with class-specific variances takes
  # one conditional-maximisation cycle per call -- so the M-step is repeated on
  # the fixed weights until the BCH-weighted log-likelihood settles. That
  # per-case log-likelihood is what the fit reports, and its variance is the
  # case-clustered sandwich over the joint parameter vector (R/step3_joint.R).
  if (inherits(model_state$sm, "nested")) {
    wll <- function(sm) rowSums(bch_resp * log_likelihood(sm, Y))
    model_state$sm <- init_params(model_state$sm, Y, bch_resp)
    prev <- -Inf
    for (iter in seq_len(max_iter)) {
      model_state$sm <- m_step(model_state$sm, Y, bch_resp, weights = weights)
      cur <- sum(weights * wll(model_state$sm))
      if (abs(cur - prev) < 1e-10 * max(1, abs(prev))) break
      prev <- cur
    }
    model_state$lower_bound <- wll(model_state$sm)
    return(.attach_nested_step3_vcov(
      model_state, Y, bch_resp, NULL, weights, se = se,
      strata = model_state$strata, cluster = model_state$cluster,
      X = X, assignment = assignment, bch_assign = A))
  }

  model_state$sm <- init_params(model_state$sm, Y, bch_resp)
  model_state$sm <- m_step(model_state$sm, Y, bch_resp)

  # Compute and store the full sandwich variance-covariance matrix of the class
  # means. This only applies to distal_continuous (which stores $means);
  # categorical models do not store $means and are excluded by the guard below.
  # Raw BCH weights are used here (before sample-weight multiplication), as
  # required by the Wald test for equality of means.
  if (inherits(model_state$sm, "distal_continuous") &&
      !is.null(model_state$sm$parameters$means)) {
    Y_vec    <- as.numeric(Y[, 1])
    K        <- model_state$n_components
    mu       <- as.vector(model_state$sm$parameters$means)

    # Restrict to cases with an observed outcome. The m_step above masks the
    # same rows, so the means this covariance matrix belongs to are estimated
    # from these cases only; Nk must be recomputed on the mask rather than over
    # all rows, or the denominators would not match the weights in the sums.
    obs     <- is.finite(Y_vec)
    bch_obs <- bch_resp[obs, , drop = FALSE]
    Y_obs   <- Y_vec[obs]
    Nk      <- colSums(bch_obs)

    if (sum(obs) < 2L || any(abs(Nk) < 1e-8)) {
      # Leave Sigma_mu unset: .wald_omnibus_means() falls back to the SEs when
      # it is NULL, whereas a matrix of zeros would be taken at face value.
      model_state$sm$parameters$Sigma_mu <- NULL
    } else {
      Sigma_mu <- matrix(0, K, K)
      for (j in 1:K) {
        for (k in j:K) {
          cov_jk <- sum(bch_obs[, j] * bch_obs[, k] *
                          (Y_obs - mu[j]) * (Y_obs - mu[k])) /
            (Nk[j] * Nk[k])
          Sigma_mu[j, k] <- cov_jk
          Sigma_mu[k, j] <- cov_jk
        }
      }
      model_state$sm$parameters$Sigma_mu <- Sigma_mu
    }
  }

  # distal_continuous_pooled's cov_theta used to be rebuilt here, and that copy
  # of the design was built over every row rather than over the rows with an
  # observed outcome that m_step actually fits -- so on data with a partly
  # missing outcome the covariance carried far more information than the
  # estimates did. m_step() now builds it alongside the coefficients, from the
  # same U, W and mask, and there is no second copy to drift.

  # A block the joint code can pack takes the joint BCH variance as its
  # one-block case: the case-clustered sandwich, plus under `se = "corrected"`
  # step 1's uncertainty carried through the rebuilt BCH weights
  # (R/step3_joint.R). A class regression or a categorical outcome takes the
  # sandwich under every `se`: the Hessian of the weighted M-step counts each
  # case fully in every class and runs 18-32% narrow. A continuous outcome
  # keeps the sandwich of its means built above unless `se = "corrected"`.
  if (.nested_packable(model_state$sm) &&
      (identical(se, "corrected") ||
       inherits(model_state$sm, c("covariate", "distal_pooled"))))
    model_state <- .attach_nested_step3_vcov(
      model_state, Y, bch_resp, NULL, weights, se = se,
      strata = model_state$strata, cluster = model_state$cluster,
      X = X, assignment = assignment, bch_assign = A)

  return(model_state)
}

# Apply the Maximum Likelihood (ML) Correction (Vermunt, 2010; Bolck et al., 2004)
#
# The core identity (Vermunt 2010, eq. 13) is:
#   P(a_i | z_i) = sum_k P(x=k | z_i) * P(a_i | x=k)
#
# P(a_i | x=k) is fixed from the proportional classification table built from
# step-1 posteriors:
#   C_prop[j,k]     = sum_i w_i * resp1[i,j] * A1[i,k]      (K x K)
#   C_row_norm[j,k] = C_prop[j,k] / sum_k C_prop[j,k]        = P(a=k | x=j)
#
# A1 is the assigned-class variable: resp1 itself under proportional assignment
# (where C_prop is symmetric), or the modal indicator under modal.
#
# For individual i under proportional (soft) assignment:
#   P(a_i | x=j) = sum_k resp1[i,k] * C_row_norm[j,k]
#
# The EM for the expanded dataset (K records per person, weight resp1[i,k]) gives:
#   Z_mat[i,k] = sum_j P(x=j | z_i) * C_row_norm[j,k]       (normaliser)
#   W[i,j]     = P(x=j | z_i) * (R %*% t(C_row_norm))[i,j]  (posterior weights)
#   where R[i,k] = resp1[i,k] / Z_mat[i,k]
#
# W sums to 1 per row and is passed to m_step as classification weights.
# Convergence LL = sum_i w_i * sum_k resp1[i,k] * log Z_mat[i,k].
fit_ml <- function(model_state, X, Y, max_iter = 1000, abs_tol = 1e-10,
                   rel_tol = 1e-10, se = "corrected",
                   assignment = "proportional") {

  if (inherits(model_state$sm, c("distal_continuous", "distal_continuous_regression"))) {
    warning(paste(
      "ML correction is not recommended for continuous distal outcomes.",
      "Bakk & Vermunt (2016) show BCH is preferred in this case,",
      "as ML requires strong distributional assumptions.",
      "Consider correction = 'BCH' instead."
    ))
  }

  # Retain every case. Missing covariates are completed under the class-invariant
  # Gaussian marginal (endogenous-constrained-x; Sterba, 2014) instead of being
  # listwise deleted, which sacrifices efficiency and biases estimates under MAR
  # (Little, 1992; Little & Zhang, 2011). A missing distal outcome is left as NA
  # and masked by the structural model's own FIML likelihood.
  Y       <- complete_structural_covariates(model_state$sm, Y)
  keep    <- rep(TRUE, nrow(Y))
  X_clean <- X[keep, , drop = FALSE]
  Y_clean <- Y[keep, , drop = FALSE]
  w_clean <- model_state$sample_weights[keep]

  # Keep the sub-model's survey design row-aligned with the retained cases so
  # that variance code running inside m_step (e.g. the covariate M-step) sees
  # strata and cluster vectors that match Y_clean.
  if (!is.null(model_state$sm$strata)) {
    model_state$sm$strata  <- model_state$sm$strata[keep]
    model_state$sm$cluster <- model_state$sm$cluster[keep]
  }

  if (nrow(X_clean) == 0) stop("All rows have missing covariates.")

  K <- model_state$n_components

  # Step 1: compute and freeze posteriors from the measurement model.
  e_res_step1 <- e_step(model_state, X_clean, NULL)
  resp_step1  <- exp(e_res_step1$log_resp)        # n_clean x K

  # The assigned-class variable whose classification error the correction
  # inverts. Under proportional assignment it is the posterior itself; under
  # modal it is the indicator of the most likely class (Bolck, Croon &
  # Hagenaars, 2004; Vermunt, 2010). Everything downstream of the frozen
  # posteriors reads A1.
  A1 <- if (assignment == "modal") get_modal_resp(resp_step1) else resp_step1

  # Classification error matrix. Rows index the true class, which is never
  # observed and so is always weighted by the posteriors; columns index the
  # assigned variable A1, which is what the rule changes. As in fit_bch(), A1
  # must be crossed with `resp_step1` and not with itself -- a modal A1 has one
  # 1 per row, so t(A1) %*% A1 is diagonal and C_row_norm would come out as the
  # identity, leaving the correction with nothing to invert. Under proportional
  # assignment A1 is `resp_step1` and this is the same matrix as before.
  #
  # ORIENTATION. `C_row_norm` is the TRANSPOSE of .classification_error()'s
  # `D`, which is the same table written with the assigned class down the rows
  # and columns summing to one. Neither is symmetric. Whichever way round a
  # caller wants it, transpose on purpose rather than by assumption.
  C_prop     <- t(resp_step1 * w_clean) %*% A1            # K x K
  Nk         <- colSums(resp_step1 * w_clean)
  C_row_norm <- sweep(C_prop, 1, Nk, "/")                 # K x K, row-normalised

  # Initialise structural model.
  model_state$sm <- init_params(model_state$sm, Y_clean, A1)
  model_state$sm <- m_step(model_state$sm, Y_clean, A1, weights = w_clean)

  # A categorical outcome is likewise a model of P(o | k), not of P(k), so its
  # class proportions are step 3's own parameters too: they start at the
  # assigned shares and are re-estimated at every M-step (Vermunt, 2010,
  # eq. 14), and the variance block below treats them as estimated.
  cat_out <- inherits(model_state$sm, c("distal_pooled", "distal_regression"))
  if (cat_out) pi_k_clean <- colSums(A1 * w_clean) / sum(w_clean)
  # A lone continuous outcome is a density f(y | k), not a model of P(k), so
  # step 3 has to carry the class proportions itself; they are re-estimated
  # with the outcome parameters at every M-step (Vermunt, 2010, eq. 14).
  cont_out <- inherits(model_state$sm, c("distal_continuous",
                                          "distal_continuous_pooled",
                                          "distal_continuous_regression"))
  if (cont_out) pi_step3 <- colSums(A1 * w_clean) / sum(w_clean)

  prev_ll <- -Inf

  step3_e <- function(sm) {

    log_sm <- log_likelihood(sm, Y_clean)

    # E-step.
    #
    # For discrete outcomes (distal_pooled, distal_regression), the correct
    # ML EM formulation following Vermunt (2010, eq. 14) uses raw outcome
    # probabilities P(o_i|x=j) = exp(log_sm)[i,j] without row-normalisation:
    #   Z_mat[i,k] = sum_j pi_j * P(o_i|x=j) * C_row_norm[j,k]
    #   W_eff[i,j] = pi_j * P(o_i|x=j) * sum_k [resp1[i,k]/Z[i,k]] * C_norm[j,k]
    #   LL         = sum_i sum_k resp1[i,k] * log Z_mat[i,k]
    #
    # Continuous outcomes (distal_continuous*) use the same formulation with
    # the proportions pi_step3 estimated in step 3. Row-normalising the
    # densities instead would be the posterior under equal class sizes, which
    # pulls every ambiguous case toward the small classes.

    if (cont_out) {
      lsh   <- .row_max(log_sm)
      joint <- sweep(exp(log_sm - lsh), 2, pi_step3, "*")
      Z_mat <- joint %*% C_row_norm
      ll_case    <- rowSums(A1 * log(pmax(Z_mat, 1e-300))) + as.vector(lsh)
      current_ll <- sum(w_clean * ll_case)
      W <- joint * ((A1 / pmax(Z_mat, 1e-300)) %*% t(C_row_norm))
    } else if (inherits(sm, c("distal_pooled", "distal_regression"))) {
      po_given_x <- exp(log_sm)                      # n_clean x K
      Z_mat      <- sweep(po_given_x, 2, pi_k_clean, "*") %*% C_row_norm
      ll_case    <- rowSums(A1 * log(pmax(Z_mat, 1e-300)))
      current_ll <- sum(w_clean * ll_case)
      RC <- A1 / pmax(Z_mat, 1e-300)
      W  <- sweep(po_given_x, 2, pi_k_clean, "*") *
        (RC %*% t(C_row_norm))
    } else {
      lsh     <- .row_max(log_sm)
      sm_prob <- exp(log_sm - lsh)
      sm_prob <- sm_prob / rowSums(sm_prob)
      Z_mat   <- sm_prob %*% C_row_norm
      ll_case    <- rowSums(A1 * log(pmax(Z_mat, 1e-300)))
      # A nested model's log_sm is log P(k | z) plus the distal log-densities,
      # which do not sum to one over k, so the row normalisation above drops
      # log sum_k of them from the log-likelihood (W is unaffected: the
      # normaliser cancels in the ratio below). Put it back so the reported
      # value is the step-3 log-likelihood itself. rowSums(A1) is 1 under
      # either assignment rule.
      if (inherits(sm, "nested"))
        ll_case <- ll_case + as.vector(lsh) +
          log(rowSums(exp(log_sm - as.vector(lsh))))
      current_ll <- sum(w_clean * ll_case)
      R <- A1 / pmax(Z_mat, 1e-300)
      W <- sm_prob * (R %*% t(C_row_norm))
    }
    list(W = W, ll_case = ll_case, current_ll = current_ll)
  }

  for (iter in seq_len(max_iter)) {

    e <- step3_e(model_state$sm)
    W <- e$W; ll_case <- e$ll_case; current_ll <- e$current_ll

    if (iter > 1) {
      change <- current_ll - prev_ll
      denom  <- max(abs(prev_ll), 1e-9)
      if (is.na(change) || abs(change) < abs_tol || abs(change / denom) < rel_tol) break
    }
    prev_ll <- current_ll

    model_state$sm <- m_step(model_state$sm, Y_clean, W, weights = w_clean)
    if (cont_out) pi_step3 <- colSums(W * w_clean) / sum(w_clean)
    if (cat_out)  pi_k_clean <- colSums(W * w_clean) / sum(w_clean)
  }

  # EM stops a lone class regression short of the maximum (gradient ~1e-2); a
  # few Newton steps on the same penalised objective finish it, so the variance
  # attached below is evaluated at the estimate it describes. The E-step is
  # re-taken so the posteriors and log-likelihood written back match.
  if (max_iter > 0 && inherits(model_state$sm, "covariate")) {
    model_state$sm$parameters$beta <- .step3_newton_polish(
      model_state$sm$parameters$beta,
      .covariate_design(model_state$sm, Y_clean), A1, C_row_norm, w_clean,
      .bayes_alpha(model_state$sm, "latent"))
    e <- step3_e(model_state$sm)
    W <- e$W; ll_case <- e$ll_case; current_ll <- e$current_ll
  }

  # With max_iter = 0 the loop never ran: take one E-step at the start so the
  # posteriors and per-case log-likelihood written back below exist.
  if (max_iter == 0) {
    e <- step3_e(model_state$sm)
    W <- e$W; ll_case <- e$ll_case
  }

  if (max_iter > 0 && inherits(model_state$sm, "distal_continuous"))
    model_state$sm <- .attach_step3_distal_vcov(
      model_state$sm, Y_clean, A1, C_row_norm, w_clean, pi_step3, se)

  # ============================================================================
  # Variance estimation for structural models under ML step-3.
  #
  # The Q-function Hessian stored by m_step reflects only the expected
  # complete-data curvature, not the marginal LL curvature, which leads to
  # underestimation of variance.  The robust sandwich estimator is used instead
  # (Bakk, Oberski & Vermunt, 2014):
  #
  #   V = B^{-1} M B^{-1}
  #
  # where B = -H_marg is the Hessian of the marginal log-likelihood
  #   L(theta) = sum_i sum_k resp1[i,k] * log Z_mat[i,k]
  # and M = sum_i s_i s_i^T is the outer product of person-level scores.
  #
  # The covariate case is delegated to R/step3_variance.R, which computes B and
  # the scores in closed form rather than by numerical differentiation and adds
  # the second variance component the two distal branches below still omit: the
  # uncertainty carried over from step 1, which is estimated and not known.
  #
  # For distal_pooled the parameters form one joint vector and the corrected
  # variance is stored in $hessian (singular), which is where downstream
  # inference reads from.
  #
  # For distal_regression the parameters are class-specific; the marginal LL
  # Hessian is block-diagonal across classes.  The per-class sandwich is
  # computed independently for each class k and stored in $hessians[[k]]
  # (plural list), which is where downstream inference reads from.
  # ============================================================================

  # --- covariate (class-membership regression) --------------------------------
  #
  # Downstream inference reads $V_robust in preference to $hessian, so attaching
  # it here is enough.
  if (inherits(model_state$sm, "covariate")) {
    model_state <- .attach_step3_covariate_vcov(
      model_state, X_clean, Y_clean, A1, C_row_norm, w_clean, se = se,
      assignment = assignment)
  }

  # A covariate block inside a nested structural model shares its step-three
  # likelihood with the distal blocks, so the single-block formulas above do not
  # apply; the variance is taken over the joint parameter vector instead
  # (R/step3_joint.R). A block that file cannot pack keeps the naive Hessian,
  # labelled as such rather than reported as though it were corrected.
  if (inherits(model_state$sm, "nested")) {
    model_state <- .attach_nested_step3_vcov(
      model_state, Y_clean, A1, C_row_norm, w_clean, se,
      strata = model_state$strata[keep], cluster = model_state$cluster[keep],
      X = X_clean, assignment = assignment)
    if (is.null(model_state$sm$parameters$vcov_joint) &&
        !is.null(model_state$sm$models[["predictor"]]) &&
        inherits(model_state$sm$models$predictor, "covariate"))
      model_state$sm$models$predictor$parameters$V_method <-
        "Q-function Hessian (uncorrected; covariate combined with a distal outcome)"
  }

  # --- distal_pooled / distal_regression (lone categorical outcome) ----------
  #
  # One sandwich over the joint vector (class-size logits, outcome
  # parameters): the sizes are estimated in step 3 too, so the outcome
  # parameters' variance is their block of it. distal_pooled stores that block
  # as one shared $hessian; distal_regression stores each class's own block in
  # $hessians[[k]], which is where downstream inference reads from.

  if (cat_out) {
    meat_fn <- function(S) .score_meat(S, model_state$strata[keep],
                                       model_state$cluster[keep],
                                       isTRUE(model_state$has_survey_design))
    if (inherits(model_state$sm, "distal_pooled")) {
      bp_dim <- dim(model_state$sm$parameters$beta_pooled)
      V <- .step3_cat_vcov(
        model_state$sm, Y_clean, A1, C_row_norm, w_clean, pi_k_clean,
        as.vector(model_state$sm$parameters$beta_pooled),
        function(sm, th) {
          sm$parameters$beta_pooled <- matrix(th, bp_dim[1], bp_dim[2]); sm },
        meat_fn)
      model_state$sm$parameters$hessian <- -pinv(V)
    } else {
      b_dim <- dim(model_state$sm$parameters$betas)          # c(K, M-1, D)
      V <- .step3_cat_vcov(
        model_state$sm, Y_clean, A1, C_row_norm, w_clean, pi_k_clean,
        as.vector(model_state$sm$parameters$betas),
        function(sm, th) { sm$parameters$betas <- array(th, b_dim); sm },
        meat_fn)
      cls <- slice.index(model_state$sm$parameters$betas, 1L)
      for (kk in seq_len(b_dim[1])) {
        ik <- which(cls == kk)
        model_state$sm$parameters$hessians[[kk]] <- -pinv(V[ik, ik, drop = FALSE])
      }
    }
  }

  # The written-back posteriors are W itself, not a fresh recombination.
  #
  # Vermunt (2010, eq. 13) treats step 3 as a genuine LC model fit on an
  # expanded data set (K pseudo-records per case, the assigned-class variable
  # A1 as a fixed-response indicator, the covariates as predictors of x); see
  # this file's header comment. W is that model's own E-step posterior
  # P(x | a_i, z_i), aggregated back to one row per case by construction --
  # it is exactly what m_step() has been fitting the structural model against
  # every iteration, at whatever iterate the loop above converged to. There is
  # no second, independent posterior to compute here.
  #
  # `keep` is always all-TRUE ("Retain every case", above), so W already
  # covers every row of X and Y -- recomputing a posterior "for the full
  # data" by re-running the E-step on the raw indicators and recombining it
  # with the fitted covariate model, as this block used to do, is not a
  # second way of getting the same answer: the LC likelihood the EM loop
  # maximises is P(a_i | z_i), not a product of the indicators' own
  # measurement-model density and the classification-error term, so that
  # recombination does not correspond to any quantity Vermunt's derivation
  # defines. `internal/validation-tests/test-step3-prior-validation.R`
  # pins colMeans(W) on its benchmark (0.4419 / 0.2954 / 0.2627, to four
  # decimals); the old recombination did not reach those figures.
  model_state$log_resp    <- log(pmax(W, 1e-300))
  model_state$lower_bound <- ll_case

  return(model_state)
}

# The ML step-3 covariance of a lone continuous outcome's class means. Step 3
# estimates the class sizes, the means and the variances together, so the
# means' variance comes from the curvature of that joint likelihood,
#   l_i = sum_s a_is log sum_k pi_k f(y_i | k) C[k, s],
# and not from the M-step's weighted means, which hold every case's class
# weights fixed and so understate it. `se = "hessian"` gives the inverse
# observed information; otherwise the sandwich H^-1 B H^-1 (White, 1982).
# Step 1's estimates are treated as known.
.attach_step3_distal_vcov <- function(sm, Y, A1, Cn, w, pk, se) {
  K   <- sm$n_components
  pk  <- pmax(pk, 1e-12)
  par <- c(log(pk[-K] / pk[K]), .qn_pack_sm(sm)$par)
  ia  <- seq_len(K - 1L)
  case_ll <- function(th) {
    a  <- c(th[ia], 0)
    p  <- exp(a - max(a)); p <- p / sum(p)
    ls <- log_likelihood(.qn_unpack_sm(sm, th[-ia]), Y)
    lsh <- .row_max(ls)
    Z  <- sweep(exp(ls - lsh), 2, p, "*") %*% Cn
    rowSums(A1 * log(pmax(Z, 1e-300))) + as.vector(lsh)
  }
  idx  <- seq_along(par)
  S    <- .qn_fd_scores(case_ll, par, idx)
  grad <- function(x) colSums(w * .qn_fd_scores(case_ll, x, idx))
  h    <- 1e-4 * pmax(1, abs(par))
  H    <- vapply(idx, function(j) {
    a <- par; a[j] <- a[j] + h[j]
    b <- par; b[j] <- b[j] - h[j]
    (grad(a) - grad(b)) / (2 * h[j])
  }, numeric(length(par)))
  Hi <- tryCatch(solve(-(H + t(H)) / 2), error = function(e) NULL)
  if (is.null(Hi)) return(sm)
  V  <- if (identical(se, "hessian")) Hi else Hi %*% crossprod(w * S) %*% Hi
  m  <- K - 1L + seq_len(K)
  sm$parameters$Sigma_mu <- V[m, m, drop = FALSE]
  sm$parameters$ses      <- matrix(sqrt(pmax(diag(V)[m], 0)), ncol = 1)
  sm$parameters$V_method <- if (identical(se, "hessian"))
    "Observed information (step 3)" else "Sandwich (step 3)"
  sm
}

# The ML step-3 sandwich of a lone categorical outcome. The per-case
# log-likelihood is l_i = sum_s a_is log sum_k pi_k P(o_i | k) C[k, s], and
# the class sizes pi are step 3's parameters as well as the outcome's, so the
# variance of the outcome parameters is their block of the sandwich over
# (size logits, outcome parameters); fixing the sizes at their estimate would
# treat them as known. `put(sm, theta)` writes an outcome vector back onto the
# structural model. Returns that block.
.step3_cat_vcov <- function(sm, Y, A1, Cn, w, pk, theta0, put, meat) {
  K   <- length(pk)
  pk  <- pmax(pk, 1e-12)
  par <- c(log(pk[-K] / pk[K]), theta0)
  ia  <- seq_len(K - 1L)
  case_ll <- function(th) {
    a <- c(th[ia], 0)
    p <- exp(a - max(a)); p <- p / sum(p)
    Z <- sweep(exp(log_likelihood(put(sm, th[-ia]), Y)), 2, p, "*") %*% Cn
    rowSums(A1 * log(pmax(Z, 1e-300)))
  }
  n_p <- length(par)
  eps <- 1e-4
  S <- matrix(0, nrow(Y), n_p)
  for (j in seq_len(n_p)) {
    e <- replace(numeric(n_p), j, eps)
    S[, j] <- (case_ll(par + e) - case_ll(par - e)) / (2 * eps)
  }
  tot <- function(th) sum(w * case_ll(th))
  H <- matrix(0, n_p, n_p)
  for (j in seq_len(n_p)) {
    for (k in j:n_p) {
      ej <- replace(numeric(n_p), j, eps)
      ek <- replace(numeric(n_p), k, eps)
      H[j, k] <- H[k, j] <-
        (tot(par + ej + ek) - tot(par + ej - ek) -
           tot(par - ej + ek) + tot(par - ej - ek)) / (4 * eps^2)
    }
  }
  B_inv <- pinv(-H)
  V <- B_inv %*% meat(S) %*% B_inv
  V[-ia, -ia, drop = FALSE]
}
