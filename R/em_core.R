# ==============================================================================
# S3 Core EM Algorithm
# ==============================================================================

# Helper to check if structural model contains a covariate
has_covariate <- function(sm) {
  if (is.null(sm)) return(FALSE)
  if (inherits(sm, "covariate")) return(TRUE)
  if (inherits(sm, "nested")) {
    return(any(sapply(sm$models, function(m) inherits(m, "covariate"))))
  }
  return(FALSE)
}

# Helper to initialize random responsibilities
initialize_resp <- function(n_samples, n_components) {
  resp <- matrix(runif(n_samples * n_components), nrow = n_samples)
  resp <- sweep(resp, 1, rowSums(resp), "/")
  return(resp)
}

# E-step: Calculate responsibilities
e_step <- function(model_state, X, Y = NULL) {
  # 1. Measurement model likelihood
  log_prob <- log_likelihood(model_state$mm, X)

  # 2. Add Structural model likelihood (if active)
  if (!is.null(Y) && !is.null(model_state$sm)) {
    log_prob <- log_prob + log_likelihood(model_state$sm, Y)
  }

  # 3. Add Marginal Prior ONLY if the covariate model is not actively providing
  # class probabilities.  A covariate SM is "active" only when Y is non-NULL AND
  # the SM is a covariate type — in that case log_likelihood(sm, Y) already
  # encodes P(class | z_i) so adding log_weights would double-count.
  # When Y = NULL (e.g. during step-1 of 3-step estimation), the covariate SM
  # contributes nothing and marginal class weights MUST still be applied.
  covariate_active <- !is.null(Y) && has_covariate(model_state$sm)
  if (!covariate_active) {
    log_weights <- log(model_state$weights + 1e-15)
    log_prob    <- sweep(log_prob, 2, log_weights, "+")
  }

  # 4. Log-sum-exp trick for stability
  log_prob_norm <- apply(log_prob, 1, function(x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  })

  log_resp <- sweep(log_prob, 1, log_prob_norm, "-")
  return(list(log_resp = log_resp, log_prob_norm = log_prob_norm))
}

# M-step: Update model parameters
m_step_core <- function(model_state, X, Y, log_resp, alpha = NULL) {
  resp <- exp(log_resp)

  # --- BAYESIAN PRIOR (Class Weights) ---
  # An explicit `alpha` still wins — fit_lta() passes its `smoothing` argument
  # down this path — but with none supplied the strength comes from the model's
  # own `bayes_constants`, which is also where the emission M-steps read theirs.
  K <- model_state$n_components
  alpha <- alpha %||% .bayes_alpha(model_state, "latent")
  prior_obs <- alpha / K

  # Sampling / frequency weights must enter the M-step, not only the reported
  # log-likelihood: a case carrying weight w contributes w times to every
  # sufficient statistic. refine_lbfgs() already weighted correctly, so on the
  # default binary/Gaussian path this only moves the starting point it is handed
  # and the fitted values barely shift. It matters where that refinement does not
  # run: polytomous and mixed measurement models, which are outside its
  # whitelist, and refine = FALSE, which is the path BLRT replicates take.
  # Passed through only when the weights are non-trivial, so unweighted fits keep
  # their exact previous numerical behaviour.
  w <- model_state$sample_weights
  weighted <- !is.null(w) && length(w) == nrow(resp) && any(w != 1)

  nk <- if (weighted) colSums(sweep(resp, 1, w, "*")) else colSums(resp)
  nk_prior <- nk + prior_obs
  model_state$weights <- nk_prior / sum(nk_prior)

  if (weighted) {
    model_state$mm <- m_step(model_state$mm, X, resp, weights = w)
  } else {
    model_state$mm <- m_step(model_state$mm, X, resp)
  }

  if (!is.null(Y) && !is.null(model_state$sm)) {
    if (weighted) {
      model_state$sm <- m_step(model_state$sm, Y, resp, weights = w)
    } else {
      model_state$sm <- m_step(model_state$sm, Y, resp)
    }
  }

  return(model_state)
}

# The emissions refine_lbfgs() knows how to polish, kept in one place because
# two separate decisions read it: whether to run the refinement at all, and how
# tightly EM itself must converge. Those two were previously decided
# independently, which is how the second one came to be wrong for every emission
# outside this list.
.refine_supported <- c("bernoulli", "bernoulli_nan",
                       "gaussian_diag", "gaussian_diag_nan",
                       "gaussian_unit", "gaussian_unit_nan")

# Would refine_lbfgs() climb from wherever EM stops to the penalised optimum for
# this emission? A block model is refinable when its sub-model is (the flat view
# maps the blocks onto one wide parameter matrix); a nested model never is,
# because its sub-models are heterogeneous and there is no single packing.
.is_refinable <- function(mm) {
  if (inherits(mm, "nested")) return(FALSE)
  view <- .refine_time_block_view(mm) %||% list(mm = mm)
  # An across-class equality constraint on the variances has no expression in
  # the refinement's parameterisation, which packs one log(sd) per class-item
  # cell; polishing such a model would walk straight off the constraint surface.
  # `col_map` ties columns together but says nothing about rows.
  if (isTRUE(view$mm$variances_equal)) return(FALSE)
  class(view$mm)[1] %in% .refine_supported
}

# The EM stopping rule for an emission that L-BFGS will not polish.
#
# The package default stops once the log-likelihood moves less than a thousandth
# of *itself*, which is deliberately loose: for the emissions the refinement
# covers, EM only has to reach the right neighbourhood and L-BFGS does the rest.
# Applied to an emission with no second stage it stops after a handful of
# iterations and that is the answer the user gets. Measured against a rule tight
# enough to be at the optimum, on simulated data with n_init = 5:
#
#   count LCA, K=3, n=800, 6 items        6 iterations,  7.0 of log-likelihood short
#   count LCA, K=4, n=2000, 8 items       6 iterations, 27.5 short
#   polytomous LCA, K=3, n=800            7 iterations,  5.1 short
#   mixed binary/continuous/count        --             8.1 short
#
# and those are not harmless decimals: a Poisson rate came back 2.8 off its
# converged value and a response probability 0.56 off, which is a different
# class profile, not a rounding difference.
#
# The value is chosen from the accuracy/cost curve rather than by taking the
# tightest rule available. At abs = 1e-4 all four models above land within 0.006
# of their maximum; tightening to 1e-8 buys the remaining 0.006 for three to
# four times the iterations. The relative term is kept only as a safety valve so
# that a very large sample cannot iterate indefinitely; at 1e-8 it does not bind
# until the log-likelihood is in the tens of thousands, and even there it is
# four orders of magnitude tighter than the default.
.em_tol_unpolished <- list(abs = 1e-4, rel = 1e-8)

# L-BFGS refinement after EM convergence — Penalised Maximum Likelihood (PM).
#
# EM converges to a Q-function fixed point, not necessarily the PM optimum.
# L-BFGS steps on the full penalised log-posterior climb past that fixed point.
#
# The objective uses a penalized maximum likelihood (PM) formulation:
#   log P(ϑ) = log L(X; ϑ) + log p(ϑ)
# where the priors (Bayes constants, all defaulting to 1 — see
# R/bayes_constants.R) are:
#   Weights  : (α_latent/K)    · Σ_k  log π_k
#   Bernoulli: (α_categorical/K) · Σ_k Σ_j [ π̂_j·log π_kj + (1-π̂_j)·log(1-π_kj) ]
#   Variances: −(α_variances/2K) · Σ_k Σ_j [ log σ²_kj + s²_j / σ²_kj ]
# with π̂_j the observed marginal probability of item j and s²_j its observed
# marginal variance ("conservative null model" in both cases).
#
# The variance term is not optional bookkeeping: it is what makes this objective
# *bounded*. Without it, L-BFGS re-optimises log(sd) with no constraint, and a
# class that has latched onto a few near-identical cases can drive one variance
# towards zero and the likelihood towards +Inf — walking straight through the
# M-step's regularisation, which is why that regularisation is now a matching
# prior rather than a floor (see m_step.gaussian_diag()). The term here and the
# one in the M-step read the same stored constant so the two stages cannot
# disagree about what is being maximised.
#
# Using the prior instead of box constraints is the correct approach: it lets
# the gradient pull parameters freely while the prior provides soft penalisation
# proportional to the evidence — which heavily stabilizes the estimation.
#
# Supported: bernoulli, bernoulli_nan, gaussian_diag, gaussian_diag_nan,
#            gaussian_unit, gaussian_unit_nan.
# No-op for nested models and other types. The missing-data (_nan) variants
# share the complete-data objective and gradient; observed-data masking is
# applied wherever the indicator matrix contributes, so refinement behaves
# consistently whether or not the data contain missing cells.
#
# Parameterisation (unconstrained):
#   bernoulli     : logit(pis) [K×J]  + log-ratio weights [K-1]
#   gaussian_diag : means [K×J] + log(sd) [K×J] + log-ratio weights [K-1]
#   gaussian_unit : means [K×J] + log-ratio weights [K-1]
#
# `.objective_only` is an internal hook, not part of the fitting path: it returns
# the packed starting values together with the value and gradient closures
# instead of running the optimiser. It exists so the analytical gradient can be
# checked against a finite-difference one directly, on the real objective rather
# than on a re-derivation of it — a wrong sign or a missing prior term still
# produces a fit that converges, so the only way to catch one is to differentiate
# the function the optimiser actually sees.
refine_lbfgs <- function(model_state, X, Y = NULL, max_iter = 500,
                         .objective_only = FALSE) {
  # A time-block (repeated-measures) model is refined on the same footing as a
  # flat one: its per-occasion parameter blocks are viewed as one wide matrix,
  # and `col_map` records which stacked column draws on which free column. Items
  # held equal across occasions therefore share a single free column, so L-BFGS
  # optimises the constrained parameterisation directly. For every other model
  # `col_map` is the identity and this reduces to the original code path.
  tb   <- .refine_time_block_view(model_state$mm)
  view <- tb %||% list(mm = model_state$mm, col_map = NULL)

  mm_type <- class(view$mm)[1]
  if (!mm_type %in% .refine_supported) return(model_state)
  if (isTRUE(view$mm$variances_equal)) return(model_state)   # see .is_refinable()
  if (inherits(view$mm, "nested")) return(model_state)
  # K=1 has no weight parameters; the M-step already gives the exact analytic
  # solution (item marginals), so L-BFGS is a no-op and the K-2 index arithmetic
  # below produces an out-of-bounds sequence that triggers a sweep() warning.
  if (model_state$n_components == 1L) return(model_state)

  # Collapse the missing-data variants onto their complete-data family so the
  # packing, likelihood, and gradient branches treat them identically. Where the
  # data contain NA, the branches below mask the missing cells (na.rm); with
  # complete data they fall through to the faster matrix-multiply forms.
  fam    <- sub("_nan$", "", mm_type)
  has_na <- anyNA(X)

  K  <- model_state$n_components
  J  <- ncol(X)
  sw <- model_state$sample_weights

  col_map <- view$col_map %||% seq_len(J)
  if (length(col_map) != J) return(model_state)
  P <- max(col_map)                       # number of free parameter columns
  tied <- P < J
  # Index sets used to fold the full-column gradient onto the free columns.
  # factor() with explicit levels keeps the groups in free-column order; a bare
  # split() on integers would order them as strings (1, 10, 2, ...).
  map_groups <- if (tied)
    split(seq_len(J), factor(col_map, levels = seq_len(P))) else NULL

  # Sum the columns of a K x J gradient within each tie group -> K x P.
  fold_cols <- function(g)
    matrix(vapply(map_groups, function(gi) rowSums(g[, gi, drop = FALSE]),
                  numeric(K)), nrow = K, ncol = P)

  # Observed marginal probabilities — the "conservative null model"
  # as the centre of the Dirichlet prior for binary items. Under an equality
  # constraint the tied columns share one prior, centred on their pooled
  # marginal, matching the prior the constrained M-step applies.
  #
  # The marginal is weighted, as in m_step.bernoulli. An unweighted mean would
  # centre the prior on a different point from the one EM used, so the same data
  # supplied as a response-pattern table with frequency weights and as expanded
  # case-level rows would not give the same answer.
  X0     <- replace(X, is.na(X), 0)
  obs_w  <- (!is.na(X)) * sw
  marginal <- colSums(X0 * sw) / pmax(colSums(obs_w), 1e-12)
  if (tied)
    marginal <- vapply(map_groups, function(g) mean(marginal[g]), numeric(1))
  marginal <- pmax(pmin(marginal, 1 - 1e-7), 1e-7)

  # Observed marginal *variance* per free column — the centre of the prior on
  # the Gaussian class variances, and the exact counterpart of `marginal` above.
  #
  # Under a tie group this pools the way the constrained M-step pools: it is the
  # variance of the *stacked* data around the *stacked* mean, not the mean of
  # the per-column variances. The two coincide only when the tied columns share
  # a mean, and where they differ, taking the average would give the polish a
  # different prior from the one EM applied — which is the same class of bug
  # this whole term exists to fix, one level down. (m_step.blocks() estimates an
  # invariant item by stacking the blocks and calling the sub-model's M-step
  # once, so the marginal it computes is the stacked one by construction.)
  alpha_var <- .bayes_alpha(view$mm, "variances")
  s2_free <- if (fam == "gaussian_diag") {
    sum_w   <- colSums(obs_w)
    sum_wx  <- colSums(X0 * sw)
    sum_wx2 <- colSums(X0^2 * sw)
    pool <- function(cols) {
      swt <- sum(sum_w[cols])
      if (!is.finite(swt) || swt <= 0) return(1)
      mu <- sum(sum_wx[cols]) / swt
      v  <- sum(sum_wx2[cols]) / swt - mu^2
      if (!is.finite(v) || v <= 0) 1 else v
    }
    if (tied) vapply(map_groups, pool, numeric(1))
    else      vapply(seq_len(J), pool, numeric(1))
  } else NULL

  # ── Pack initial parameters ────────────────────────────────────────────────
  wts <- pmax(model_state$weights, 1e-15)
  log_ratio_w <- log(wts[-K] / wts[K])   # K-1 free weight params (last anchored = 0)

  # Free-column view of the starting values (identical to the full matrix when
  # nothing is tied; one representative column per group otherwise).
  free_cols <- if (tied) vapply(map_groups, function(g) g[1L], integer(1)) else
    seq_len(J)

  if (fam == "bernoulli") {
    pis  <- pmax(pmin(view$mm$parameters$pis[, free_cols, drop = FALSE],
                      1 - 1e-7), 1e-7)
    par0 <- c(qlogis(as.vector(pis)), log_ratio_w)   # K*P + K-1

  } else if (fam == "gaussian_diag") {
    means <- as.vector(view$mm$parameters$means[, free_cols, drop = FALSE])
    sds   <- sqrt(as.vector(view$mm$parameters$covariances[, free_cols, drop = FALSE]))
    par0  <- c(means, log(pmax(sds, 1e-7)), log_ratio_w)   # 2*K*P + K-1

  } else {  # gaussian_unit
    par0 <- c(as.vector(view$mm$parameters$means[, free_cols, drop = FALSE]),
              log_ratio_w)                                  # K*P + K-1
  }

  n_obs <- sum(sw)   # effective sample size (supports survey weights)

  # ── Penalised log-posterior + analytical gradient ──────────────────────────
  # Providing an analytical gradient to optim() eliminates the O(p) function
  # evaluations per iteration that numerical finite-difference gradient
  # estimation requires (p = K*J + K-1 ≈ 47 for K=4, J=11).  Without it,
  # L-BFGS spends ~50 * (p+1) ≈ 5700 function calls per fit; with it, the
  # inner loop reduces to ~50 * 4 ≈ 200 calls — a ~25x speedup that makes
  # refine_lbfgs fast at any sample size.
  #
  # The gradient derivations follow from the chain rule on the PM objective
  # PM(θ) = log L(X;θ) + log p(θ) and the reparameterisations:
  #   · Bernoulli pis on logit scale:  ∂PM/∂logit(π_kj) = Σ_i sw_i R_ik (x_ij − π_kj) + (m̂_j − π_kj)/K
  #   · Class weights on log-ratio:   ∂PM/∂lr_k = Σ_i sw_i (R_ik − w_k) + 1/K − w_k  (k < K)
  #   · Gaussian means (unit var):     ∂PM/∂μ_kj = Σ_i sw_i R_ik (x_ij − μ_kj)
  #   · Gaussian means (diag):         ∂PM/∂μ_kj = Σ_i sw_i R_ik (x_ij − μ_kj) / σ_kj²
  #   · Gaussian log-sd (diag):        ∂PM/∂log(σ_kj) = Σ_i sw_i R_ik [(x_ij−μ_kj)²/σ_kj² − 1]

  # Helper: given current params, compute neg-PM value AND gradient in one pass.
  neg_pm_and_grad <- function(par) {

    # ── Decode class weights ─────────────────────────────────────────────────
    n_w  <- K - 1L
    lr   <- c(par[(length(par) - n_w + 1):length(par)], 0)
    lr   <- lr - (max(lr) + log(sum(exp(lr - max(lr)))))
    log_w <- lr
    w_vec <- exp(log_w)

    # ── Decode measurement model ─────────────────────────────────────────────
    # Parameters are carried on the free columns and expanded to all stacked
    # columns through col_map, which is the identity unless items are tied.
    expand <- function(m) if (tied) m[, col_map, drop = FALSE] else m

    if (fam == "bernoulli") {
      pis_free <- pmax(pmin(matrix(plogis(par[seq_len(K * P)]), K, P),
                            1 - 1e-15), 1e-15)
      pis_p    <- expand(pis_free)
      if (has_na) {
        # Per-component sum with missing cells dropped (FIML). Slower than the
        # matrix multiply but the only NA-safe form.
        log_lik <- matrix(0, nrow(X), K)
        for (k in seq_len(K))
          log_lik[, k] <- rowSums(
            sweep(X, 2, log(pis_p[k, ]), "*") +
              sweep(1 - X, 2, log(1 - pis_p[k, ]), "*"),
            na.rm = TRUE)
      } else {
        # Vectorised log-likelihood: n×K  (matrix multiply — no R for-loop over k)
        log_lik <- X %*% t(log(pis_p)) + (1 - X) %*% t(log(1 - pis_p))
      }

    } else if (fam == "gaussian_diag") {
      means_free <- matrix(par[seq_len(K * P)],             K, P)
      sds_free   <- pmax(matrix(exp(par[(K*P+1):(2*K*P)]),  K, P), 1e-7)
      means_p <- expand(means_free)
      sds_p   <- expand(sds_free)
      log_lik <- matrix(0, nrow(X), K)
      for (k in seq_len(K))
        log_lik[, k] <- rowSums(
          dnorm(X, mean = matrix(means_p[k,], nrow(X), J, byrow = TRUE),
                sd   = matrix(sds_p[k,],  nrow(X), J, byrow = TRUE), log = TRUE),
          na.rm = TRUE)

    } else {  # gaussian_unit
      means_free <- matrix(par[seq_len(K * P)], K, P)
      means_p    <- expand(means_free)
      log_lik <- matrix(0, nrow(X), K)
      for (k in seq_len(K))
        log_lik[, k] <- rowSums(
          dnorm(X, mean = matrix(means_p[k,], nrow(X), J, byrow = TRUE),
                sd = 1, log = TRUE), na.rm = TRUE)
    }

    # ── Posterior responsibilities (needed for both value and gradient) ───────
    log_joint <- sweep(log_lik, 2, log_w, "+")
    mx        <- apply(log_joint, 1, max)
    log_norm  <- mx + log(rowSums(exp(sweep(log_joint, 1, mx, "-"))))
    R         <- exp(sweep(log_joint, 1, log_norm, "-"))   # n × K posteriors

    # ── Observed-data log-likelihood (normalised by n_obs for scale-invariance) ─
    obs_ll <- sum(sw * log_norm) / n_obs

    # ── Priors (also normalised by n_obs) ─────────────────────────────────────
    log_prior_w <- (1 / K) * sum(log_w) / n_obs
    log_prior_pis <- if (fam == "bernoulli")
      sum((marginal / K) %*% t(log(pis_free)) +
            ((1 - marginal) / K) %*% t(log(1 - pis_free))) / n_obs
    else 0

    # Truncated inverse-Wishart on the class variances (see the header). Written
    # on the free columns, so a tied item contributes its prior once — matching
    # the constrained M-step, which estimates it once from the stacked blocks.
    log_prior_var <- if (fam == "gaussian_diag" && alpha_var > 0) {
      var_free <- sds_free^2
      -0.5 * (alpha_var / K) *
        sum(log(var_free) + sweep(1 / var_free, 2, s2_free, "*")) / n_obs
    } else 0

    val <- -(obs_ll + log_prior_w + log_prior_pis + log_prior_var)

    # ── Analytical gradient ───────────────────────────────────────────────────
    swR <- sweep(R, 1, sw, "*")          # sw_i * R_ik,  n × K
    nk  <- colSums(swR)                  # Σ_i sw_i R_ik, length K

    grad <- numeric(length(par))

    # The data part of each gradient is computed on all stacked columns and then
    # folded onto the free columns; tied columns therefore accumulate the score
    # from every occasion they cover. The prior is added once per free column.
    if (fam == "bernoulli") {
      if (has_na) {
        # Missing cells inform neither item j's score nor its effective count,
        # so zero the missing contributions and count only observed cells per
        # item. With complete data nk_obs[k, ] == nk[k] and this reduces to the
        # matrix-multiply form below.
        obs    <- !is.na(X)
        X0     <- X; X0[!obs] <- 0
        nk_obs <- t(swR) %*% obs
        g_data <- t(swR) %*% X0 - pis_p * nk_obs
      } else {
        g_data <- t(swR) %*% X - sweep(pis_p, 1, nk, "*")
      }
      if (tied) g_data <- fold_cols(g_data)
      g_pis <- g_data + (matrix(marginal, K, P, byrow = TRUE) - pis_free) / K
      grad[seq_len(K * P)] <- as.vector(-g_pis) / n_obs

    } else if (fam == "gaussian_diag") {
      g_mu  <- matrix(0, K, J)
      g_lsd <- matrix(0, K, J)
      for (k in seq_len(K)) {
        res_k <- sweep(X, 2, means_p[k,], "-")
        z2    <- sweep(res_k^2, 2, sds_p[k,]^(-2), "*")          # (x − μ)² / σ²
        g_mu[k,]  <- colSums(swR[,k] * sweep(res_k, 2, sds_p[k,]^(-2), "*"),
                             na.rm = has_na)
        g_lsd[k,] <- colSums(swR[,k] * (z2 - 1), na.rm = has_na)
      }
      if (tied) { g_mu <- fold_cols(g_mu); g_lsd <- fold_cols(g_lsd) }
      # Prior score, added once per free column (the data score above is summed
      # over every stacked column the free column covers; the prior is not).
      #   ∂/∂log σ_kj of −(α/2K)[log σ²_kj + s²_j/σ²_kj] = (α/K)(s²_j/σ²_kj − 1)
      if (alpha_var > 0)
        g_lsd <- g_lsd + (alpha_var / K) *
          (sweep(1 / sds_free^2, 2, s2_free, "*") - 1)
      grad[seq_len(K * P)]  <- as.vector(-g_mu)  / n_obs
      grad[(K*P+1):(2*K*P)] <- as.vector(-g_lsd) / n_obs

    } else {  # gaussian_unit
      g_mu <- matrix(0, K, J)
      for (k in seq_len(K))
        g_mu[k,] <- colSums(swR[,k] * sweep(X, 2, means_p[k,], "-"), na.rm = has_na)
      if (tied) g_mu <- fold_cols(g_mu)
      grad[seq_len(K * P)] <- as.vector(-g_mu) / n_obs
    }

    g_w <- nk[seq_len(n_w)] - n_obs * w_vec[seq_len(n_w)] +
      1/K - w_vec[seq_len(n_w)]
    grad[(length(par) - n_w + 1):length(par)] <- -g_w / n_obs

    list(val = val, grad = grad)
  }

  # Wrapper returning just value (for optim fn=)
  neg_pm_val  <- function(par) neg_pm_and_grad(par)$val
  # Wrapper returning just gradient (for optim gr=)
  neg_pm_grad <- function(par) neg_pm_and_grad(par)$grad

  if (isTRUE(.objective_only))
    return(list(par0 = par0, fn = neg_pm_val, gr = neg_pm_grad, fam = fam,
                K = K, P = P, tied = tied))

  # ── Run L-BFGS with analytical gradient (unconstrained) ────────────────────
  fit <- tryCatch(
    optim(par0, neg_pm_val, gr = neg_pm_grad, method = "L-BFGS-B",
          control = list(maxit = max_iter, factr = 1e7)),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence > 1) return(model_state)

  # ── Unpack ────────────────────────────────────────────────────────────────
  par <- fit$par
  n_w <- K - 1L

  lr   <- c(par[(length(par) - n_w + 1):length(par)], 0)
  lr   <- lr - (max(lr) + log(sum(exp(lr - max(lr)))))
  model_state$weights <- exp(lr)

  expand_out <- function(m) if (tied) m[, col_map, drop = FALSE] else m
  refined <- list()

  if (fam == "bernoulli") {
    refined$pis <- expand_out(
      pmax(pmin(matrix(plogis(par[1:(K*P)]), nrow = K, ncol = P), 1-1e-7), 1e-7))

  } else if (fam == "gaussian_diag") {
    refined$means       <- expand_out(matrix(par[1:(K*P)], nrow = K, ncol = P))
    refined$covariances <- expand_out(
      matrix(exp(par[(K*P+1):(2*K*P)]), nrow = K, ncol = P)^2)

  } else {  # gaussian_unit
    refined$means <- expand_out(matrix(par[1:(K*P)], nrow = K, ncol = P))
  }

  if (is.null(tb)) {
    for (nm in names(refined)) model_state$mm$parameters[[nm]] <- refined[[nm]]
  } else {
    # Scatter the wide refined matrix back into the per-occasion sub-models.
    model_state$mm <- .refine_time_block_write(model_state$mm, refined)
  }

  # Re-run E-step with the refined parameters so log_resp and lower_bound
  # reflect the true post-refinement posteriors.
  e_res <- e_step(model_state, X, Y)
  model_state$log_resp    <- e_res$log_resp
  model_state$lower_bound <- e_res$log_prob_norm

  return(model_state)
}

# Run the EM loop for a single random initialization.
#
# `init_state` resumes from a state some earlier call already produced instead
# of drawing fresh starting values, which is what makes the two-stage search in
# fit_em() possible: a short first pass ranks the starts, and only the survivors
# are run on to convergence from exactly where they stopped.
fit_single_init <- function(model_state, X, Y, max_iter = 1000,
                            abs_tol = 1e-3, rel_tol = 1e-3, refine = TRUE,
                            init_state = NULL) {
  n_samples <- nrow(X)

  # The default stopping rule is deliberately loose because the emissions it was
  # written for are polished afterwards by refine_lbfgs(), which climbs from
  # wherever EM stopped to the penalised optimum. An emission outside that
  # whitelist has no such second stage, so where EM stops is the answer, and it
  # gets the tighter rule automatically. An emission may still name its own rule
  # by carrying `em_tol`, which takes precedence over both.
  em_tol <- model_state$mm$em_tol
  if (is.null(em_tol) && !.is_refinable(model_state$mm))
    em_tol <- .em_tol_unpolished
  if (!is.null(em_tol)) {
    abs_tol <- em_tol$abs
    rel_tol <- em_tol$rel
  }

  if (is.null(init_state)) {
    model_state$weights <- rep(1 / model_state$n_components,
                               model_state$n_components)

    # Initialize parameters directly
    model_state$mm <- init_params(model_state$mm, X, NULL)
    if (!is.null(Y) && !is.null(model_state$sm)) {
      model_state$sm <- init_params(model_state$sm, Y, NULL)
    }
  } else {
    model_state$weights <- init_state$weights
    model_state$mm      <- init_state$mm
    if (!is.null(init_state$sm)) {
      model_state$sm <- init_state$sm
    } else if (!is.null(Y) && !is.null(model_state$sm)) {
      # A resumed state produced by the two-stage search always carries its own
      # structural model, but a warm start built from a measurement-only fit
      # (see `warm_start` in fit_em()) has nothing to say about it, and an
      # emission that never saw init_params() has no parameters for m_step_core
      # to update.
      model_state$sm <- init_params(model_state$sm, Y, NULL)
    }
  }

  # Initialize scalar trackers for convergence
  prev_total_ll <- -Inf
  converged <- FALSE
  n_iter <- 0

  for (iter in 1:max_iter) {
    # 1. E-STEP
    e_res <- e_step(model_state, X, Y)

    # Store the FULL VECTOR in the model state (for weights/BIC)
    log_prob_vector <- e_res$log_prob_norm

    # Calculate the SCALAR TOTAL LL for the convergence check
    # This uses the weights we added to support survey data!
    current_total_ll <- sum(model_state$sample_weights * log_prob_vector)

    # 2. CONVERGENCE CHECK (Using scalars)
    if (iter > 1) {
      change <- current_total_ll - prev_total_ll

      if (!is.na(change) && (abs(change) < abs_tol || abs(change / max(abs(prev_total_ll), 1e-9)) < rel_tol)) {
        converged <- TRUE
        n_iter <- iter
        model_state$lower_bound <- log_prob_vector # Final vector storage
        break
      }
    }

    # 3. M-STEP: Update parameters
    model_state <- m_step_core(model_state, X, Y, e_res$log_resp)

    # Prepare for next iteration
    prev_total_ll <- current_total_ll
    n_iter <- iter
  }

  model_state$converged <- converged
  model_state$n_iter <- n_iter
  model_state$lower_bound <- log_prob_vector # Ensure the vector is returned
  model_state$log_resp <- e_res$log_resp

  # L-BFGS refinement per restart so all restarts compete on PM likelihood.
  # With analytical gradients this is fast (~0.3s at n=5000), so the n_init-fold
  # cost is acceptable and selection bias from EM-only ranking is avoided.
  # Skipped when refine = FALSE (BLRT bootstrap replicates).
  if (isTRUE(refine)) model_state <- refine_lbfgs(model_state, X, Y)

  return(model_state)
}

# Multi-start EM
#
# An emission whose EM is slow may ask for a two-stage search by carrying
# `em_stage1`, and for a higher iteration ceiling by carrying `em_max_iter`.
# Both are opt-in: an emission that declares neither runs exactly the loop it
# always did.
#
# `warm_start` is a function of (model_state, X, Y) returning a state to run as
# one extra restart alongside the random ones, or NULL to add none. Random
# starts search the whole parameter space uniformly, which is the wrong prior
# for a model whose parameters are split into per-group blocks that only mean
# the same thing when their classes are aligned; a start built from a fitted
# restriction of the model is in the aligned basin by construction. It competes
# on log-likelihood like any other restart, so a warm start that turns out to be
# a poor one costs a restart and changes nothing else.
fit_em <- function(model_state, X, Y, n_init = 1, max_iter = 1000,
                   random_state = NULL, refine = TRUE, warm_start = NULL) {

  # An emission may raise the ceiling on itself, but never lower one the caller
  # asked for: `max_iter` is a documented argument of fit_mixture().
  em_max <- model_state$mm$em_max_iter
  if (!is.null(em_max)) max_iter <- max(max_iter, em_max)

  run_from <- function(state)
    fit_single_init(model_state, X, Y, max_iter = max_iter, refine = refine,
                    init_state = state)

  # The warm start gets the tight (absolute) rule even where the emission is one
  # L-BFGS polishes. The default rule is *relative to the log-likelihood's own
  # magnitude*, so in the thousands it fires at a change of several units, and
  # the warm start has several *hundred* units of climbing to do along a ridge —
  # measured on the multi-country validation model it stopped after 6 iterations,
  # 7 units short, and the polish then converged to the nearest local optimum
  # rather than following the ridge. A random start does not have this problem
  # because it begins near nothing in particular.
  #
  # This is not a quirk of one model. Biernacki, Celeux and Govaert (2003,
  # p. 568), setting up the standard comparison of EM initialisation strategies:
  # "We do not use stopping criteria based on the relative change of the
  # estimates or loglikelihood because the slow convergence of the EM makes such
  # criteria hazardous." Their own short runs stop on progress relative to
  # progress *made so far*, (L^q - L^{q-1}) / (L^q - L^0) — scale-free in the way
  # that matters, where a rule relative to |L| is not. An absolute rule is the
  # cheaper fix and is what `.em_tol_unpolished` already provides.
  #
  # An emission that names its own `em_tol` still overrides this inside
  # fit_single_init().
  run_warm <- function(state)
    fit_single_init(model_state, X, Y, max_iter = max_iter, refine = refine,
                    init_state = state,
                    abs_tol = .em_tol_unpolished$abs,
                    rel_tol = .em_tol_unpolished$rel)

  ll_of <- function(s) sum(s$sample_weights * s$lower_bound)

  # Built once, before either search below, so that a failure to build one is
  # reported (and tolerated) in a single place. The warm start is deterministic,
  # so it is unaffected by the seed handling in the loops.
  warm <- NULL
  if (is.function(warm_start)) warm <- warm_start(model_state, X, Y)

  stage <- model_state$mm$em_stage1
  if (!is.null(stage) && n_init > 1L) {
    # Stage 1 — rank the starting values cheaply.
    #
    # Where each iteration is expensive and hundreds are needed, running every
    # restart to convergence spends nearly all of its time climbing hills that
    # are then discarded. A short pass separates the promising basins from the
    # hopeless ones for a fraction of the cost, and only the survivors are run
    # on. This is the standard two-stage multi-start scheme, used for exactly the
    # same reason, and the survivors resume from where stage 1 left them rather
    # than restarting, so nothing is thrown away.
    candidates <- vector("list", n_init)
    for (init in seq_len(n_init)) {
      if (!is.null(random_state)) set.seed(random_state + init)
      candidates[[init]] <- fit_single_init(model_state, X, Y,
                                            max_iter = stage$iter,
                                            refine = FALSE)
    }
    lls  <- vapply(candidates, ll_of, numeric(1))
    keep <- order(lls, decreasing = TRUE)[
      seq_len(min(n_init, max(stage$min_keep, ceiling(stage$frac * n_init))))]

    best_model <- NULL
    best_total_ll <- -Inf
    for (i in keep) {
      # The seed is restored so that the second stage of a given start draws the
      # same random numbers it would have drawn in a single-stage run; nothing
      # below uses them today, but a stochastic M-step later would.
      if (!is.null(random_state)) set.seed(random_state + i)
      fitted_state <- run_from(candidates[[i]])
      current_ll   <- ll_of(fitted_state)
      if (is.null(best_model) || current_ll > best_total_ll) {
        best_total_ll <- current_ll
        best_model    <- fitted_state
      }
    }
    # The warm start skips stage 1: ranking it against random starts on a short
    # pass would defeat its purpose, since the basin it seeds is the one a short
    # pass is worst at recognising (the GMM measurements in the roadmap).
    if (!is.null(warm)) {
      fitted_state <- run_warm(warm)
      if (ll_of(fitted_state) > best_total_ll) best_model <- fitted_state
    }
    return(best_model)
  }

  best_model <- NULL
  best_total_ll <- -Inf

  for (init in 1:n_init) {
    if (!is.null(random_state)) set.seed(random_state + init)

    # Each restart includes L-BFGS refinement so selection is by PM likelihood.
    fitted_state <- fit_single_init(model_state, X, Y, max_iter = max_iter,
                                    refine = refine)

    current_ll <- ll_of(fitted_state)
    if (is.null(best_model) || current_ll > best_total_ll) {
      best_total_ll <- current_ll
      best_model    <- fitted_state
    }
  }

  if (!is.null(warm)) {
    fitted_state <- run_warm(warm)
    if (is.null(best_model) || ll_of(fitted_state) > best_total_ll)
      best_model <- fitted_state
  }

  return(best_model)
}
