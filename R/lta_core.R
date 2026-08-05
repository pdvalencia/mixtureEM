# ==============================================================================
# Latent Transition Analysis - EM engine (forward-backward / Baum-Welch)
# ==============================================================================
#
# In a latent transition model a person occupies a latent *status* at each
# occasion and may move between statuses. Writing S_1..S_T for the statuses,
#
#   P(y_i) = Σ_{s_1..s_T} δ_{s_1} Π_{t=2..T} τ_{t, s_t | s_{t-1}}
#                                 Π_{t=1..T} Π_j ρ_{j t | s_t}(y_ijt)
#
# with δ the latent status prevalences at Time 1, τ_t the transition matrix from
# occasion t-1 to t, and ρ the item-response parameters (Collins & Lanza, 2010,
# sec. 7.5).
#
# The sum runs over K^T status sequences and is never evaluated directly: the
# E-step is the forward-backward recursion, and EM built on it is the Baum-Welch
# algorithm. This model cannot use the engine in em_core.R, whose e_step()
# assumes a single K-valued latent variable with an n x K posterior. Here the
# E-step yields a posterior per occasion (γ) and a posterior over consecutive
# status pairs (ξ); it is ξ that identifies the transitions.
#
# The measurement models are reused rather than rebuilt. The emission term at
# occasion t is log_likelihood(<measurement model>, X_t), so binary, polytomous,
# continuous and mixed indicators all work unchanged, and the "_nan" variants
# supply FIML: an occasion at which a case has no observed item contributes a
# flat term, which is the correct treatment of wave attrition. Across-time
# equality constraints on ρ use the same time_blocks code as fit_rmlca().

# ------------------------------------------------------------------------------
# E-step
# ------------------------------------------------------------------------------

# Emission log-densities, one n x K matrix per occasion.
.lta_emission_loglik <- function(mm, X) {
  J <- mm$n_items
  lapply(seq_len(mm$n_times), function(t)
    log_likelihood(mm$models[[t]], X[, .time_block_cols(t, J), drop = FALSE]))
}

# Forward-backward in log space.
#
# Returns the per-case log-likelihood, the occasion-wise posteriors γ, and the
# weighted pairwise transition counts Ξ_t[k, l] = Σ_i w_i P(S_t = k, S_{t+1} = l | y_i).
# The pairwise term is accumulated one (k, l) pair at a time: each element is a
# probability, so its logarithm is non-positive and exponentiating is safe,
# while never materialising an n x K x K array.
#
# `log_delta` is either a length-K vector or, when covariates predict the initial
# status, an n x K matrix. Each element of `log_tau` is either a K x K matrix or,
# when covariates predict transitions, a list of K matrices of size n x K (one
# per origin status). The three helpers below hide that distinction so the
# recursions are written once.
.lta_forward_backward <- function(logB, log_delta, log_tau, weights,
                                  keep_pairwise = FALSE) {
  Tn <- length(logB)
  n  <- nrow(logB[[1]])
  K  <- ncol(logB[[1]])

  # Log-transition into destination l, as an n x K matrix over origin statuses.
  tau_to <- function(t, l) {
    LT <- log_tau[[t]]
    if (is.list(LT))
      do.call(cbind, lapply(seq_len(K), function(k) LT[[k]][, l]))
    else matrix(LT[, l], n, K, byrow = TRUE)
  }
  # Log-transitions out of origin k, as an n x K matrix over destinations.
  tau_from <- function(t, k) {
    LT <- log_tau[[t]]
    if (is.list(LT)) LT[[k]] else matrix(LT[k, ], n, K, byrow = TRUE)
  }

  la <- vector("list", Tn)
  lb <- vector("list", Tn)

  la[[1]] <- if (is.matrix(log_delta)) logB[[1]] + log_delta else
    sweep(logB[[1]], 2, log_delta, "+")
  if (Tn > 1L) for (t in 2:Tn) {
    prev <- la[[t - 1]]
    m <- matrix(0, n, K)
    for (l in seq_len(K))
      m[, l] <- logsumexp(prev + tau_to(t - 1L, l), MARGIN = 1)
    la[[t]] <- m + logB[[t]]
  }

  lb[[Tn]] <- matrix(0, n, K)
  if (Tn > 1L) for (t in (Tn - 1):1) {
    nxt <- logB[[t + 1]] + lb[[t + 1]]
    m <- matrix(0, n, K)
    for (k in seq_len(K))
      m[, k] <- logsumexp(nxt + tau_from(t, k), MARGIN = 1)
    lb[[t]] <- m
  }

  ll <- logsumexp(la[[Tn]], MARGIN = 1)          # length n

  gamma <- lapply(seq_len(Tn), function(t)
    exp(la[[t]] + lb[[t]] - ll))

  xi   <- vector("list", max(Tn - 1L, 0L))
  pair <- if (keep_pairwise) vector("list", max(Tn - 1L, 0L)) else NULL
  if (Tn > 1L) for (t in seq_len(Tn - 1L)) {
    emit_next <- logB[[t + 1]] + lb[[t + 1]]     # n x K
    Xi_t <- matrix(0, K, K)
    pair_t <- if (keep_pairwise) vector("list", K) else NULL
    for (k in seq_len(K)) {
      base_k <- la[[t]][, k] - ll
      block  <- exp(sweep(emit_next, 1, base_k, "+") + tau_from(t, k))
      Xi_t[k, ] <- colSums(block * weights)
      if (keep_pairwise) pair_t[[k]] <- block
    }
    xi[[t]] <- Xi_t
    if (keep_pairwise) pair[[t]] <- pair_t
  }

  list(ll = ll, gamma = gamma, xi = xi, pairwise = pair)
}

# ------------------------------------------------------------------------------
# M-step pieces
# ------------------------------------------------------------------------------

# Normalise a vector of expected counts into probabilities with a Dirichlet
# prior of total mass `alpha` spread over the admissible cells, then re-impose
# any structural zeros. LTA on sparse tables reaches boundary solutions
# routinely; this is the standard smoothing remedy (Collins & Lanza, sec. 6.11).
.lta_normalise <- function(counts, alpha, allowed = NULL) {
  if (is.null(allowed)) allowed <- rep(TRUE, length(counts))
  if (!any(allowed)) return(rep(1 / length(counts), length(counts)))
  out <- numeric(length(counts))
  prior <- alpha / sum(allowed)
  num <- counts[allowed] + prior
  out[allowed] <- num / sum(num)
  out
}

# ------------------------------------------------------------------------------
# Latent classes above the chain
# ------------------------------------------------------------------------------
#
# A mixture latent Markov model puts a second latent variable above the chain: C
# classes, each with its own initial distribution and its own transition
# matrices, sharing one measurement model (Vermunt, "Mover-Stayer Models"; van
# de Pol & Langeheine, 1990),
#
#   P(y_i) = Σ_c π_c · P(y_i | class c)
#
# where each P(y_i | class c) is the ordinary latent Markov likelihood already
# implemented above. The mover-stayer model is the case where one class's
# transition matrix is the identity, which needs no new machinery at all: the
# `tau_allowed` mask that carries `forbidden_transitions` already forces it,
# since .lta_normalise() over a single admissible cell returns exactly 1.
#
# The E-step is therefore *unchanged* - .lta_forward_backward() is simply run
# once per class. `.lta_class_state()` below is what makes that literal: it
# swaps in one class's parameters and hands back something with the shape the
# single-chain code already expects, so .lta_log_delta() and .lta_log_tau() do
# not know that classes exist.
#
# The per-class parameters live in `delta_c`, `tau_c` and `tau_allowed_c`, which
# are lists indexed by class. `.lta_pack()` writes them back out to `delta`,
# `tau` and `tau_allowed` at the end of a fit: for a single class those are the
# vector and list of matrices they have always been, so nothing downstream
# changes, and for several they gain a leading class index.

.lta_class_state <- function(state, c) {
  state$delta       <- state$delta_c[[c]]
  state$tau         <- state$tau_c[[c]]
  state$tau_allowed <- state$tau_allowed_c[[c]]
  state
}

.lta_pack <- function(state) {
  if (state$n_classes == 1L) {
    state$delta       <- state$delta_c[[1]]
    state$tau         <- state$tau_c[[1]]
    state$tau_allowed <- state$tau_allowed_c[[1]]
    state$class_weights <- NULL
  } else {
    state$delta <- do.call(rbind, state$delta_c)
    dimnames(state$delta) <- list(
      paste("Class", seq_len(state$n_classes)),
      paste("Status", seq_len(state$n_statuses)))
    state$tau         <- state$tau_c
    state$tau_allowed <- state$tau_allowed_c
    names(state$class_weights) <- paste("Class", seq_len(state$n_classes))
  }
  state
}

# ------------------------------------------------------------------------------
# One EM run from the current parameter values
# ------------------------------------------------------------------------------
.lta_em <- function(state, X, max_iter = 1000, tol = 1e-8, alpha = 1.0) {
  w  <- state$weights_vec
  Tn <- state$n_times
  K  <- state$n_statuses
  C  <- state$n_classes
  n  <- nrow(X)

  prev_ll <- -Inf
  converged <- FALSE
  n_iter <- 0L

  has_delta_cov <- !is.null(state$Z_delta)
  has_tau_cov   <- !is.null(state$Z_tau)
  # With several classes the transition counts have to be re-aggregated with the
  # class posterior folded into the weights, so the recursion is asked to keep
  # its unweighted pairwise blocks rather than the totals it would form itself.
  keep_pair <- has_tau_cov || C > 1L

  # One E-step across all classes: the per-class posteriors, the class
  # posterior, and the case log-likelihood of the mixture.
  e_step <- function(state) {
    logB <- .lta_emission_loglik(state$mm, X)
    es <- lapply(seq_len(C), function(c) {
      sub <- .lta_class_state(state, c)
      .lta_forward_backward(logB, .lta_log_delta(sub), .lta_log_tau(sub), w,
                            keep_pairwise = keep_pair)
    })
    if (C == 1L)
      return(list(es = es, post = matrix(1, n, 1L), ll = es[[1]]$ll))

    lp <- vapply(seq_len(C), function(c) es[[c]]$ll, numeric(n))
    lp <- sweep(lp, 2, log(pmax(state$class_weights, 1e-300)), "+")
    ll <- logsumexp(lp, MARGIN = 1)
    list(es = es, post = exp(lp - ll), ll = ll)
  }

  for (iter in seq_len(max_iter)) {
    E <- e_step(state)
    cur_ll <- sum(w * E$ll)

    if (iter > 1L) {
      change <- cur_ll - prev_ll
      if (!is.na(change) &&
          (abs(change) < tol ||
           abs(change / max(abs(prev_ll), 1e-9)) < tol)) {
        converged <- TRUE
        n_iter <- iter
        break
      }
    }
    prev_ll <- cur_ll
    n_iter  <- iter

    # --- class weights --------------------------------------------------------
    if (C > 1L)
      state$class_weights <- .lta_normalise(colSums(E$post * w), alpha)

    for (c in seq_len(C)) {
      es <- E$es[[c]]
      # Everything this class learns is weighted by the posterior probability
      # that a case belongs to it; with one class that factor is 1 throughout
      # and the arithmetic below is the single-chain M-step unchanged.
      wc <- if (C == 1L) w else w * E$post[, c]

      # --- initial status prevalences ----------------------------------------
      # The covariate M-steps set `delta` and `tau` themselves - they are the
      # average of a case-level regression rather than a normalised count - so
      # their answers are mirrored into the per-class slots that everything
      # downstream now reads. Covariates and several classes are mutually
      # exclusive (fit_lta() refuses the combination), so `c` is always 1 here.
      if (has_delta_cov) {
        state <- .lta_mstep_delta_cov(state, es$gamma[[1]])
        state$delta_c[[c]] <- state$delta
      } else {
        state$delta_c[[c]] <- .lta_normalise(colSums(es$gamma[[1]] * wc), alpha)
      }

      # --- transition matrices ------------------------------------------------
      if (Tn > 1L && has_tau_cov) {
        state <- .lta_mstep_tau_cov(state, es)
        state$tau_c[[c]] <- state$tau
      } else if (Tn > 1L) {
        Xi <- .lta_pair_counts(es, wc, K, Tn, C)
        allowed <- state$tau_allowed_c[[c]]
        if (isTRUE(state$tau_homogeneous)) {
          pooled <- Reduce(`+`, Xi)
          tau1 <- matrix(0, K, K)
          for (k in seq_len(K))
            tau1[k, ] <- .lta_normalise(pooled[k, ], alpha, allowed[[1]][k, ])
          state$tau_c[[c]] <- rep(list(tau1), Tn - 1L)
        } else {
          for (t in seq_len(Tn - 1L)) {
            m <- matrix(0, K, K)
            for (k in seq_len(K))
              m[k, ] <- .lta_normalise(Xi[[t]][k, ], alpha, allowed[[t]][k, ])
            state$tau_c[[c]][[t]] <- m
          }
        }
      }
    }

    # --- measurement model ---------------------------------------------------
    # A list of responsibilities, one per occasion: the invariance constraints
    # in m_step.time_blocks pool them across occasions where ρ is held equal.
    # The measurement model is shared by the classes, so the responsibility it
    # sees is the class-mixed status posterior P(S_t = k | y_i), which is what
    # the additivity of the complete-data log-likelihood calls for.
    state$mm <- m_step(state$mm, X, .lta_mixed_gamma(E, Tn, C),
                       weights = if (all(w == 1)) NULL else w, alpha = alpha)
  }

  # Final E-step so the stored posteriors match the returned parameters.
  E <- e_step(state)

  state$ll_case   <- E$ll
  state$loglik    <- sum(w * E$ll)
  state$gamma     <- .lta_mixed_gamma(E, Tn, C)
  state$xi        <- Reduce(function(a, b) Map(`+`, a, b),
                            lapply(seq_len(C), function(c)
                              .lta_pair_counts(E$es[[c]],
                                               if (C == 1L) w else w * E$post[, c],
                                               K, Tn, C)))
  state$converged <- converged
  state$n_iter    <- n_iter

  if (C > 1L) {
    state$class_posterior <- E$post
    state$gamma_by_class  <- lapply(E$es, `[[`, "gamma")
    state$xi_by_class <- lapply(seq_len(C), function(c)
      .lta_pair_counts(E$es[[c]], w * E$post[, c], K, Tn, C))
  }
  .lta_pack(state)
}

# Weighted transition counts Ξ_t[k, l]. When the recursion was asked to keep its
# pairwise blocks they are unweighted, so the weights are applied here; when it
# was not, it has already formed the totals under the same weights.
.lta_pair_counts <- function(es, weights, K, Tn, C) {
  if (Tn < 2L) return(list())
  if (is.null(es$pairwise)) return(es$xi)
  lapply(seq_len(Tn - 1L), function(t) {
    m <- matrix(0, K, K)
    for (k in seq_len(K))
      m[k, ] <- colSums(es$pairwise[[t]][[k]] * weights)
    m
  })
}

# P(S_t = k | y_i), summing the per-class status posteriors over the class
# posterior. With one class this is that class's posterior unchanged.
.lta_mixed_gamma <- function(E, Tn, C) {
  if (C == 1L) return(E$es[[1]]$gamma)
  lapply(seq_len(Tn), function(t)
    Reduce(`+`, lapply(seq_len(C), function(c)
      E$es[[c]]$gamma[[t]] * E$post[, c])))
}

# ------------------------------------------------------------------------------
# Random starts
# ------------------------------------------------------------------------------

# Draw a probability vector from a symmetric Dirichlet, using only the gamma
# trick on rexp/rgamma-free base R: exponential variates normalised.
.rdirichlet1 <- function(K, concentration = 1) {
  g <- stats::rgamma(K, shape = concentration, rate = 1)
  if (sum(g) <= 0) return(rep(1 / K, K))
  g / sum(g)
}

.lta_random_start <- function(state, X) {
  K  <- state$n_statuses
  Tn <- state$n_times
  C  <- state$n_classes

  # Each restart begins from probabilities rather than from regression
  # coefficients; the first M-step fits the regressions to them.
  state$delta_beta <- NULL
  state$tau_beta   <- NULL

  if (C > 1L) state$class_weights <- .lta_normalise(.rdirichlet1(C, 5) * 10, 0)

  for (c in seq_len(C)) {
    state$delta_c[[c]] <- .lta_normalise(.rdirichlet1(K, 5) * 10, 0)

    for (t in seq_len(max(Tn - 1L, 0L))) {
      m <- matrix(0, K, K)
      for (k in seq_len(K)) {
        # Start with mass concentrated on staying put: transition matrices in
        # practice are diagonally dominant, and a diffuse start makes the label
        # of a status drift between occasions before the measurement model
        # settles. Across classes the concentration is varied rather than
        # shared, because classes that start alike stay alike - a mixture whose
        # chains are drawn from one distribution collapses onto a single chain.
        stay <- if (C == 1L) 2 else 2 * c / C
        draw <- .rdirichlet1(K, 1) + (seq_len(K) == k) * stay
        m[k, ] <- .lta_normalise(draw, 0, state$tau_allowed_c[[c]][[t]][k, ])
      }
      state$tau_c[[c]][[t]] <- m
    }
    if (isTRUE(state$tau_homogeneous) && Tn > 1L)
      state$tau_c[[c]] <- rep(state$tau_c[[c]][1], Tn - 1L)
  }

  state$mm <- init_params(state$mm, X, NULL)
  state
}
