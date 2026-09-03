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
# prior of total mass `alpha` spread evenly over the admissible cells, then
# re-impose any structural zeros. LTA on sparse tables reaches boundary
# solutions routinely; this is the standard smoothing remedy (Collins & Lanza,
# sec. 6.11).
#
# Two choices here are deliberate and both are sourced.
#
# The mass is one pseudo-case per conditional table, not per cell and not per
# row. Chung, Lanza & Loken (2008) put a Dirichlet with constant
# hyper-parameter 1/L on the joint class-membership probabilities - "adding just
# one observation to each class at time 1" - and under a prior of that shape the
# boundary solutions ML produces at n = 100 disappear with lower RMSE. The
# add-one-per-cell family (Jeffreys, Laplace, Goodman's add-2) is markedly worse
# on sparse transition structures, which is Fienberg & Holland's (1973) own
# finding about adding 2 to every cell of a large sparse table. Until 2026-08-31
# the mass here was one pseudo-case per *row*, which is K times what the program
# this package calibrates against puts on a transition matrix; `patterns` below
# is what fixes that, and `### 41.1 re-opened` in internal/ROADMAP.md records
# the reconciliation that measured it.
#
# The spread is even, not proportional to the destinations' marginal. That is
# the opposite of the measurement model's prior (m_step.bernoulli(), which
# centres on the item's observed marginal) and the difference is not an
# inconsistency. Fienberg & Holland (1973, sec. 6) shrink a cross-classified
# table toward its independence fit, and their risk contours show that target
# winning near an odds ratio of 1 and losing ground as the table departs from
# independence. A transition matrix departs from independence about as far as a
# table can: the destination marginal is dominated by the prevalent status, so
# shrinking a rare origin row toward it asserts that everyone moves there. Even
# is uninformative; marginal would be confidently wrong. The dilution argument
# that makes the marginal form right for the item-response probabilities does
# not apply either - a transition row has K cells however many items, occasions
# or cases the model has.
#
# `patterns` is the number of distinct predictor patterns the conditional
# distribution being normalised is estimated for -- 1 for a marginal
# distribution, K for a transition matrix's rows, K * C when a latent class sits
# above the chain as well. `alpha` is the mass for the WHOLE conditional table,
# so a row of it receives `alpha / patterns` and a cell of that row
# `alpha / (patterns * sum(allowed))`. This is the other program's constant: its
# prior on a latent variable is `alpha / (K * U0)` per cell, with U0 the count of
# unique covariate and predictor patterns, and in a chained LTA the previous
# status is a predictor. Until 2026-08-31 `patterns` was effectively 1 here, so a
# transition matrix received K times the mass it should have -- see
# `### 41.1 re-opened` in internal/ROADMAP.md, which pins the constant by
# evaluating the kernel at that program's own converged solution.
.lta_normalise <- function(counts, alpha, allowed = NULL, patterns = 1L) {
  if (is.null(allowed)) allowed <- rep(TRUE, length(counts))
  if (!any(allowed)) return(rep(1 / length(counts), length(counts)))
  out <- numeric(length(counts))
  prior <- alpha / (patterns * sum(allowed))
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
# Global (Viterbi) decoding
# ------------------------------------------------------------------------------
#
# .lta_forward_backward() answers "what is the marginal status at occasion t?",
# one occasion at a time; the most probable status at every occasion need not
# be a sequence the model gives any probability at all, if the transition
# between two locally-favoured statuses happens to be forbidden or merely
# unlikely. Viterbi decoding instead finds the single most probable *sequence*
# S_1..S_T, replacing the forward recursion's sum over origin statuses with a
# max: same tau_from() shape as .lta_forward_backward(), max instead of
# logsumexp forward, then a backward trace of which origin fed each maximum.
.lta_viterbi_one <- function(logB, log_delta, log_tau) {
  Tn <- length(logB)
  n  <- nrow(logB[[1]])
  K  <- ncol(logB[[1]])

  tau_from <- function(t, k) {
    LT <- log_tau[[t]]
    if (is.list(LT)) LT[[k]] else matrix(LT[k, ], n, K, byrow = TRUE)
  }
  row_max <- function(m) do.call(pmax, as.data.frame(m))

  lp <- vector("list", Tn)
  lp[[1]] <- if (is.matrix(log_delta)) logB[[1]] + log_delta else
    sweep(logB[[1]], 2, log_delta, "+")
  if (Tn > 1L) for (t in 2:Tn) {
    prev <- lp[[t - 1]]
    m <- matrix(0, n, K)
    for (l in seq_len(K)) {
      cand <- matrix(0, n, K)
      for (k in seq_len(K)) cand[, k] <- prev[, k] + tau_from(t - 1L, k)[, l]
      m[, l] <- row_max(cand)
    }
    lp[[t]] <- m + logB[[t]]
  }

  path <- matrix(0L, n, Tn)
  path[, Tn] <- max.col(lp[[Tn]], ties.method = "first")
  if (Tn > 1L) for (t in (Tn - 1L):1L) {
    nxt  <- path[, t + 1L]
    cand <- matrix(0, n, K)
    for (k in seq_len(K))
      cand[, k] <- lp[[t]][, k] + tau_from(t, k)[cbind(seq_len(n), nxt)]
    path[, t] <- max.col(cand, ties.method = "first")
  }
  list(path = path, logp = row_max(lp[[Tn]]))
}

# With more than one latent class, decoding has to consider class and path
# jointly rather than plugging the marginal modal class into its own chain's
# best path: the two questions are not separable, since which chain a case
# most plausibly followed depends on which class it is in. Running the single-
# chain recursion once per class and taking the per-case maximum over
# (class, path) together is exact and needs no further convention.
.lta_viterbi <- function(state) {
  C    <- state$n_classes %||% 1L
  logB <- .lta_emission_loglik(state$mm, state$data)

  if (C == 1L) {
    v <- .lta_viterbi_one(logB, .lta_log_delta(state), .lta_log_tau(state))
    return(list(path = v$path, logp = v$logp, class = NULL))
  }

  per <- lapply(seq_len(C), function(c) {
    st <- .lta_class_state(state, c)
    v  <- .lta_viterbi_one(logB, .lta_log_delta(st), .lta_log_tau(st))
    v$logp <- v$logp + log(pmax(state$class_weights[c], 1e-300))
    v
  })
  L   <- do.call(cbind, lapply(per, `[[`, "logp"))
  cls <- max.col(L, ties.method = "first")
  path <- matrix(0L, nrow(L), state$n_times)
  for (c in seq_len(C)) {
    i <- cls == c
    if (any(i)) path[i, ] <- per[[c]]$path[i, , drop = FALSE]
  }
  list(path = path, logp = L[cbind(seq_len(nrow(L)), cls)], class = cls)
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

  prev_obj <- -Inf
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
  e_step <- function(state, keep = keep_pair) {
    if (!is.null(state$ri)) return(.lta_ri_e_step(state, X, w))
    logB <- .lta_emission_loglik(state$mm, X)
    es <- lapply(seq_len(C), function(c) {
      sub <- .lta_class_state(state, c)
      .lta_forward_backward(logB, .lta_log_delta(sub), .lta_log_tau(sub), w,
                            keep_pairwise = keep)
    })
    if (C == 1L)
      return(list(es = es, post = matrix(1, n, 1L), ll = es[[1]]$ll))

    lp <- vapply(seq_len(C), function(c) es[[c]]$ll, numeric(n))
    lp <- sweep(lp, 2, log(pmax(state$class_weights, 1e-300)), "+")
    ll <- logsumexp(lp, MARGIN = 1)
    list(es = es, post = exp(lp - ll), ll = ll)
  }

  # EM is stopped on the objective it climbs -- the plain log-likelihood plus
  # the M-step's own prior -- not on the log-likelihood alone. Where that
  # objective cannot be written down exactly (.lta_log_prior() says so by
  # returning NA) the rule falls back to the plain log-likelihood. The item
  # marginals behind the measurement term are formed once here; support cannot
  # change during the loop, so it is settled once too.
  prior_marg <- if (alpha > 0 || .bayes_alpha(state$mm$models[[1]],
                                              "categorical") > 0)
    tryCatch(.lta_prior_marginals(state, X), error = function(e) NULL) else NULL
  obj_of <- function(st, ll) {
    lp <- .lta_log_prior(st, X, alpha, prior_marg)
    if (is.na(lp)) ll else ll + lp
  }
  use_pen <- !is.na(.lta_log_prior(state, X, alpha, prior_marg))

  for (iter in seq_len(max_iter)) {
    E <- e_step(state)
    cur_ll <- sum(w * E$ll)
    # `state` here is still the point that produced `E`, so the prior and the
    # log-likelihood are read at the same parameters.
    cur_obj <- if (use_pen) obj_of(state, cur_ll) else cur_ll

    if (iter > 1L) {
      change <- cur_obj - prev_obj
      if (!is.na(change) &&
          (abs(change) < tol ||
           abs(change / max(abs(prev_obj), 1e-9)) < tol)) {
        converged <- TRUE
        n_iter <- iter
        break
      }
    }
    prev_obj <- cur_obj
    n_iter  <- iter

    # --- class weights --------------------------------------------------------
    if (C > 1L)
      state$class_weights <- .lta_normalise(colSums(E$post * w), alpha)

    # A tied initial-status distribution pools the raw (pre-normalisation)
    # occasion-1 responsibilities across every class into one shared table,
    # the same "sum sufficient statistics, normalise once, broadcast back"
    # trick `tau_homogeneous` below already does across occasions instead of
    # across classes. Computed once here, ahead of the per-class loop, since
    # it needs every class's `gamma[[1]]` before any of them can be assigned.
    if (isTRUE(state$tie_initial_status) && C > 1L && !has_delta_cov) {
      raw <- Reduce(`+`, lapply(seq_len(C), function(c) {
        wc <- w * E$post[, c]
        colSums(E$es[[c]]$gamma[[1]] * wc)
      }))
      delta_tied <- .lta_normalise(raw, alpha, patterns = 1L)
    }

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
      } else if (isTRUE(state$tie_initial_status) && C > 1L) {
        state$delta_c[[c]] <- delta_tied
      } else {
        # With a latent class above the chain the initial-status distribution is
        # itself conditional on that class, so the prior mass is shared over the
        # C patterns. With one class this is the marginal case and `patterns` is
        # 1, which is what it has always been.
        state$delta_c[[c]] <- .lta_normalise(colSums(es$gamma[[1]] * wc), alpha,
                                             patterns = C)
      }

      # --- transition matrices ------------------------------------------------
      if (Tn > 1L && has_tau_cov) {
        state <- .lta_mstep_tau_cov(state, es)
        state$tau_c[[c]] <- state$tau
      } else if (Tn > 1L) {
        Xi <- .lta_pair_counts(es, wc, K, Tn, C)
        allowed <- state$tau_allowed_c[[c]]
        # A transition matrix is one conditional table over K origin patterns
        # (times C, when a latent class conditions it too), so `alpha` is its
        # mass in total rather than each row's. See .lta_normalise().
        n_pat <- K * C
        if (isTRUE(state$tau_homogeneous)) {
          pooled <- Reduce(`+`, Xi)
          tau1 <- matrix(0, K, K)
          for (k in seq_len(K))
            tau1[k, ] <- .lta_normalise(pooled[k, ], alpha, allowed[[1]][k, ],
                                        patterns = n_pat)
          state$tau_c[[c]] <- rep(list(tau1), Tn - 1L)
        } else {
          for (t in seq_len(Tn - 1L)) {
            m <- matrix(0, K, K)
            for (k in seq_len(K))
              m[k, ] <- .lta_normalise(Xi[[t]][k, ], alpha, allowed[[t]][k, ],
                                       patterns = n_pat)
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
    #
    # The measurement model's prior is NOT `smoothing`. `smoothing` is the
    # Dirichlet mass on the status and transition probabilities, which are counts
    # over K destinations; the measurement model's is a marginal-preserving prior
    # on the item-response probabilities, which is a different prior of a
    # different shape on a different table (Chung, Lanza & Loken 2008 make
    # exactly this split, sec. 3). Passing `alpha` here won the `%||%` in
    # m_step.bernoulli() and silently replaced `bayes_constants$categorical`, so
    # `?fit_lta`'s documented division of labour between the two arguments was
    # not the one in force, and `smoothing = 0` also stripped the measurement
    # prior - reintroducing the boundary estimates on ρ that the prior is there
    # to prevent.
    if (!is.null(state$ri)) {
      state <- .lta_ri_mstep(state, X, E, alpha)
    } else {
      state$mm <- m_step(state$mm, X, .lta_mixed_gamma(E, Tn, C),
                         weights = if (all(w == 1)) NULL else w)
    }
  }

  # A random intercept is identified only up to a reflection/relabelling; fix
  # it here, before the final E-step, so that step's `pq` is already computed
  # under the normalised parameterisation and never needs its own permutation.
  if (!is.null(state$ri)) state <- .lta_ri_sign_normalise(state)

  # Final E-step so the stored posteriors match the returned parameters.
  E <- e_step(state, keep = TRUE)
  if (!is.null(state$ri)) state$ri$pq <- E$ri$pq

  # Absolute entropy of the joint status-path posterior, the numerator of
  # metrics$entropy. The smoothed path posterior of a hidden Markov chain is
  # itself Markov, so this is a sum of one marginal and T-1 conditional
  # entropies over quantities the E-step has already formed -- O(n T K^2), not
  # K^T. With several classes it is the class-posterior-weighted average of
  # the per-class path entropies; the separate class entropy stays in
  # metrics$class_entropy and is deliberately not folded in.
  state$abs_ent_path <- local({
    total <- 0
    for (c in seq_len(C)) {
      wc  <- if (C == 1L) w else w * E$post[, c]
      g   <- E$es[[c]]$gamma
      h_i <- rowSums(-g[[1]] * log(g[[1]] + 1e-15))
      if (Tn > 1L) for (t in seq_len(Tn - 1L)) for (k in seq_len(K)) {
        Pk <- E$es[[c]]$pairwise[[t]][[k]]
        h_i <- h_i + rowSums(-Pk * log((Pk + 1e-15) / (g[[t]][, k] + 1e-15)))
      }
      total <- total + sum(wc * h_i)
    }
    total
  })

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

# ------------------------------------------------------------------------------
# Post-EM refinement: L-BFGS on the penalised log-posterior
# ------------------------------------------------------------------------------
#
# EM converges to a Q-function fixed point, which need not be the optimum of the
# objective itself. R/em_core.R has said so, and taken an L-BFGS step past it,
# since the mixture engine was written; the LTA driver is separate -- it needs
# the forward-backward recursion the mixture engine has no place for -- and had
# no equivalent, so it was the one estimation path in the package where EM's
# last mile was walked rather than jumped. On a near-flat ridge that mile is
# long: on the LTA-FAQ benchmark plain EM stops ~0.005 log-likelihood units
# short at the default tolerance and needs ~1,480 iterations at tol = 1e-13 to
# close the gap.
#
# The objective is the SAME penalised log-posterior the M-step maximises, not
# the plain log-likelihood. If the two stages disagreed about what is being
# maximised the polish would pull every fit off the estimator the package
# documents -- the failure R/em_core.R guards against by having its M-step and
# its refinement read one stored constant. Here the three penalties are:
#
#   delta      (alpha / K)          * sum_k log delta_k
#   tau row k  (alpha / |allowed_k|)* sum_{l allowed} log tau[k, l]
#   rho        (alpha_cat / K)      * sum_kj [ m_j log rho_kj
#                                              + (1 - m_j) log(1 - rho_kj) ]
#   log_sd     -(alpha_var / 2K)    * sum_k [ log sigma^2_kj + s^2_j / sigma^2_kj ]
#
# matching .lta_normalise() and m_step.bernoulli() term for term, with m_j the
# weighted observed marginal of item j over exactly the occasions that item's
# M-step pools (all of them for an invariant item, one for a free one -- see
# m_step.blocks()); the log_sd term matches m_step.gaussian_diag()'s own
# truncated inverse-Wishart prior the same way, with s^2_j the item's observed
# marginal variance over those occasions. Gaussian means carry no prior. With
# `smoothing = 0` and `bayes_constants$categorical = 0` (and `variances = 0`)
# every term vanishes and the objective is the plain log-likelihood.
#
# Two properties are asserted in the tests rather than argued here: the
# analytic gradient agrees with central finite differences, and at an EM fixed
# point the gradient of this objective is ~0 -- which is the proof that the
# objective is EM's own and not a near neighbour of it.

# Column layout of the unconstrained vector, in .lta_score_matrix()'s block
# order. Each entry says how to read a slice of the vector back into `state`.
#
# With several classes above the chain, a class-mixing block comes first and
# every delta/tau block gains a `c` field naming which class's own slot it
# reads (`state$delta_c[[c]]`, `state$tau_c[[c]][[i_mat]]`) -- these per-class
# lists are what `.lta_pack()` always maintains, collapsed to `state$delta`/
# `state$tau` only as a C == 1 convenience, so reading them directly here
# needs no branch on C at all: at C == 1 there is one class and `b$c` is
# always 1.
.lta_par_layout <- function(state) {
  K  <- state$n_statuses
  Tn <- state$n_times
  C  <- state$n_classes %||% 1L
  out <- list()
  if (C > 1L) out[[1L]] <- list(kind = "class", len = C - 1L)

  # Covariate models (`delta_beta`/`tau_beta`) replace the plain probability
  # tables with case-level regressions; they are mutually exclusive with
  # `C > 1` in fit_lta(), so `c` is always 1 wherever these fire.
  idx <- if (isTRUE(state$tau_homogeneous)) 1L else seq_len(Tn - 1L)
  # A tied initial-status distribution is one shared (K-1)-vector, not C of
  # them, so it gets one block for the whole model rather than one per class.
  # Covariate models replace delta with a regression and are mutually exclusive
  # with C > 1 in fit_lta(), so they never reach this.
  tied <- isTRUE(state$tie_initial_status) && C > 1L &&
    is.null(state$delta_beta)
  for (c in seq_len(C)) {
    if (!is.null(state$delta_beta)) {
      D <- ncol(state$Z_delta)
      out[[length(out) + 1L]] <- list(kind = "delta_beta", c = c,
                                      len = (K - 1L) * D)
    } else if (!tied || c == 1L) {
      out[[length(out) + 1L]] <- list(kind = "delta", c = c, len = K - 1L,
                                      tied = tied)
    }

    if (!is.null(state$tau_beta)) {
      by_origin <- identical(state$transition_effects, "by_origin")
      for (i_mat in idx) {
        if (by_origin) {
          D <- ncol(state$Z_tau)
          for (k in seq_len(K))
            out[[length(out) + 1L]] <- list(kind = "tau_beta", c = c,
                                            i_mat = i_mat, k = k,
                                            len = (K - 1L) * D)
        } else {
          D <- ncol(.lta_tau_design(state, 1L))
          out[[length(out) + 1L]] <- list(kind = "tau_beta", c = c,
                                          i_mat = i_mat, len = (K - 1L) * D)
        }
      }
    } else {
      for (i_mat in idx) for (k in seq_len(K)) {
        allowed <- which(state$tau_allowed_c[[c]][[i_mat]][k, ])
        if (length(allowed) < 2L) next
        out[[length(out) + 1L]] <- list(kind = "tau", c = c, i_mat = i_mat, k = k,
                                        allowed = allowed, len = length(allowed) - 1L)
      }
    }
  }

  J   <- state$n_items
  fam <- class(state$mm$models[[1]])[1]
  kind <- if (fam %in% c("bernoulli", "bernoulli_nan")) "rho" else "mu"
  inv <- .lta_invariant_items(state)
  for (j in seq_len(J)) {
    grps <- if (j %in% inv) list(seq_len(Tn)) else lapply(seq_len(Tn), identity)
    for (grp in grps) {
      out[[length(out) + 1L]] <- list(kind = kind, j = j, grp = grp, len = K)
      if (kind == "mu" &&
          !is.null(state$mm$models[[grp[1]]]$parameters$covariances))
        out[[length(out) + 1L]] <- list(kind = "log_sd", j = j, grp = grp, len = K)
    }
  }
  out
}

# state -> vector. Multinomial logits anchored on the last (admissible)
# category, logit for rho, identity for Gaussian means. The class-mixing
# block, where present, is anchored on the last class the same way delta is
# anchored on the last status.
.lta_par_pack <- function(state, layout) {
  K <- state$n_statuses
  C <- state$n_classes %||% 1L
  unlist(lapply(layout, function(b) {
    switch(b$kind,
      class = {
        p <- pmax(state$class_weights, 1e-12)
        log(p[seq_len(C - 1L)]) - log(p[C])
      },
      delta = {
        p <- pmax(state$delta_c[[b$c]], 1e-12)
        log(p[seq_len(K - 1L)]) - log(p[K])
      },
      delta_beta = as.vector(t(state$delta_beta[seq_len(K - 1L), , drop = FALSE])),
      tau = {
        p <- pmax(state$tau_c[[b$c]][[b$i_mat]][b$k, b$allowed], 1e-12)
        log(p[-length(p)]) - log(p[length(p)])
      },
      tau_beta = {
        B <- .lta_tau_beta(state, b$i_mat, b$k %||% 1L)
        as.vector(t(B[seq_len(K - 1L), , drop = FALSE]))
      },
      rho = {
        p <- state$mm$models[[b$grp[1]]]$parameters$pis[, b$j]
        stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12))
      },
      mu = state$mm$models[[b$grp[1]]]$parameters$means[, b$j],
      log_sd = 0.5 * log(state$mm$models[[b$grp[1]]]$parameters$covariances[, b$j]))
  }), use.names = FALSE)
}

# vector -> state. The inverse of .lta_par_pack(), writing an invariant item's
# one parameter block to every occasion that shares it.
.lta_par_unpack <- function(par, state, layout) {
  K  <- state$n_statuses
  Tn <- state$n_times
  C  <- state$n_classes %||% 1L
  pos <- 0L
  for (b in layout) {
    v <- par[pos + seq_len(b$len)]
    pos <- pos + b$len
    if (b$kind == "class") {
      p <- exp(c(v, 0) - max(c(v, 0)))
      state$class_weights <- p / sum(p)
    } else if (b$kind == "delta") {
      p <- exp(c(v, 0) - max(c(v, 0)))
      p <- p / sum(p)
      # Tied: one free vector, written to every class -- the vector->state
      # inverse of the M-step's own "normalise once, broadcast back" (the
      # `delta_tied` block in .lta_em(), R/lta_core.R).
      if (isTRUE(b$tied)) for (cc in seq_len(C)) state$delta_c[[cc]] <- p
      else state$delta_c[[b$c]] <- p
    } else if (b$kind == "tau") {
      p <- exp(c(v, 0) - max(c(v, 0)))
      p <- p / sum(p)
      row <- numeric(K)
      row[b$allowed] <- p
      state$tau_c[[b$c]][[b$i_mat]][b$k, ] <- row
    } else if (b$kind == "delta_beta") {
      D <- ncol(state$Z_delta)
      state$delta_beta <- rbind(matrix(v, K - 1L, D, byrow = TRUE), 0)
    } else if (b$kind == "tau_beta") {
      by_origin <- identical(state$transition_effects, "by_origin")
      D <- if (by_origin) ncol(state$Z_tau) else ncol(.lta_tau_design(state, 1L))
      B <- rbind(matrix(v, K - 1L, D, byrow = TRUE), 0)
      if (by_origin) state$tau_beta[[b$i_mat]][[b$k]] <- B
      else state$tau_beta[[b$i_mat]] <- B
    } else if (b$kind == "rho") {
      for (tt in b$grp)
        state$mm$models[[tt]]$parameters$pis[, b$j] <- stats::plogis(v)
    } else if (b$kind == "log_sd") {
      for (tt in b$grp)
        state$mm$models[[tt]]$parameters$covariances[, b$j] <- exp(2 * v)
    } else {
      for (tt in b$grp)
        state$mm$models[[tt]]$parameters$means[, b$j] <- v
    }
  }
  # Covariate models don't write state$delta_c/tau_c above -- they replace
  # those tables with regressions -- so the plain probability views everything
  # downstream reads (transition_matrix(), prevalence output, and delta_c/
  # tau_c themselves) are recomputed here as the case-weighted average
  # implied by the just-unpacked coefficients, the same tail
  # .lta_mstep_delta_cov()/.lta_mstep_tau_cov() run after every M-step
  # (R/lta_covariates.R:110-118, :186-195).
  if (!is.null(state$delta_beta)) {
    w <- state$weights_vec
    state$delta <- colSums(exp(.lta_log_delta(state)) * w) / sum(w)
    state$delta_c[[1L]] <- state$delta
  }
  if (!is.null(state$tau_beta)) {
    w  <- state$weights_vec
    lt <- .lta_log_tau(state)
    state$tau <- lapply(lt, function(per_k)
      t(vapply(per_k, function(M) colSums(exp(M) * w) / sum(w), numeric(K))))
    state$tau_c[[1L]] <- state$tau
  }
  # A homogeneous transition matrix is stored once per interval, so each
  # class's single free matrix has to be broadcast back over them.
  if (isTRUE(state$tau_homogeneous) && Tn > 2L)
    for (c in seq_len(C))
      state$tau_c[[c]] <- rep(state$tau_c[[c]][1], Tn - 1L)
  .lta_pack(state)
}

# Case-level log-likelihood at an arbitrary point of the packed vector. The
# plain likelihood, never the penalised objective: `smoothing` and
# `bayes_constants` enter only through the M-step's priors, and nothing on this
# path calls them. That is what the MLR scaling factor and the
# observed-information sandwich are both defined on.
#
# Deliberately not routed through .lta_score_matrix(): that function runs the
# forward-backward pass with `keep_pairwise = TRUE` and builds the whole score
# matrix, and the finite-difference Hessian that calls this does so thousands of
# times and wants neither.
.lta_ll_case <- function(state, X, par, layout) {
  st <- .lta_par_unpack(par, state, layout)
  C  <- st$n_classes %||% 1L
  logB <- .lta_emission_loglik(st$mm, X)
  if (C == 1L)
    return(.lta_forward_backward(logB, .lta_log_delta(st), .lta_log_tau(st),
                                 st$weights_vec)$ll)

  # A mixture's per-case log-likelihood is log Σ_c π_c P(y_i | class c); each
  # class runs its own chain over the shared emissions.
  lls <- vapply(seq_len(C), function(c) {
    sub <- .lta_class_state(st, c)
    .lta_forward_backward(logB, log(pmax(sub$delta, 1e-300)),
                          lapply(sub$tau, function(m) log(pmax(m, 1e-300))),
                          st$weights_vec)$ll
  }, numeric(nrow(X)))
  logsumexp(sweep(lls, 2, log(pmax(st$class_weights, 1e-300)), "+"),
            MARGIN = 1)
}

# The weighted observed marginal of item j over the occasions its M-step pools,
# which is the centre of m_step.bernoulli()'s prior for that block.
.lta_rho_prior_marginal <- function(state, X, b) {
  J <- state$n_items
  w <- state$weights_vec
  num <- 0; den <- 0
  for (tt in b$grp) {
    xj  <- X[, .time_block_cols(tt, J)[b$j]]
    obs <- !is.na(xj)
    num <- num + sum(w[obs] * xj[obs])
    den <- den + sum(w[obs])
  }
  if (den == 0) 0.5 else num / den
}

# The weighted observed marginal variance of item j over the occasions its
# M-step pools, stacking those occasions into one column and reusing
# .marginal_var() (R/gaussian.R) -- the exact quantity m_step.gaussian_diag()
# centres its own variance prior on, so the polish's prior and the M-step's
# prior are the same number by construction, not by coincidence.
.lta_var_prior_marginal <- function(state, X, b) {
  J <- state$n_items
  w <- state$weights_vec
  xs <- numeric(0); ws <- numeric(0)
  for (tt in b$grp) {
    xj <- X[, .time_block_cols(tt, J)[b$j]]
    xs <- c(xs, xj)
    ws <- c(ws, w)
  }
  .marginal_var(matrix(xs, ncol = 1), ws)
}

# The observed item marginals every measurement prior is centred on, one per
# (occasion, item), with an invariant item's pooled value repeated down its
# column. Identical to what .lta_rho_prior_marginal() returns block by block,
# but formed once: it depends only on `X` and the case weights, neither of which
# moves during EM.
.lta_prior_marginals <- function(state, X) {
  J <- state$n_items; Tn <- state$n_times; w <- state$weights_vec
  inv <- .lta_invariant_items(state)
  M <- matrix(NA_real_, Tn, J)
  for (j in seq_len(J)) {
    grps <- if (j %in% inv) list(seq_len(Tn)) else lapply(seq_len(Tn), identity)
    for (g in grps) {
      num <- 0; den <- 0
      for (tt in g) {
        xj  <- X[, .time_block_cols(tt, J)[j]]
        obs <- !is.na(xj)
        num <- num + sum(w[obs] * xj[obs])
        den <- den + sum(w[obs])
      }
      M[g, j] <- if (den == 0) 0.5 else num / den
    }
  }
  M
}

# The log-prior the M-step actually adds: the quantity whose sum with the plain
# log-likelihood EM is guaranteed not to decrease. The EM stopping rule and the
# best-restart ranking both compare that sum, because a rule that stops or ranks
# on a quantity the iteration is not climbing can fire early, late or on noise
# -- the plain log-likelihood is not monotone under a penalised M-step. See
# `### 41.1 re-opened` in internal/ROADMAP.md for the measurement, including why
# the reverse (ranking on the plain log-likelihood) looked better under a
# scoring rule that was itself the plain log-likelihood.
#
# Returns NA_real_ where the exact quantity cannot be written down, and the
# caller then falls back to the plain log-likelihood -- the rule the package
# used everywhere before 2026-08-31. An INCOMPLETE penalty would be worse than
# no penalty, because the sum would not be monotone either. The gap is the
# measurement families whose prior is not the marginal-preserving Beta one:
# Gaussian variances, Poisson rates and the polytomous Dirichlet. `smoothing`
# is inert wherever a covariate M-step replaces a normalised count, so those
# blocks are skipped on exactly the condition the M-step branches on.
# The delta/tau (and, with several classes, class-mixing) terms of the prior,
# unaffected by whether the measurement side is a plain table or a random
# intercept: factored out so `.lta_log_prior()` and `.lta_ri_log_prior()`
# (R/lta_ri.R) cannot drift into two different formulae for the same terms.
.lta_log_prior_dt <- function(state, alpha) {
  K  <- state$n_statuses
  C  <- state$n_classes %||% 1L
  Tn <- state$n_times
  val <- 0
  if (alpha <= 0) return(val)

  if (C > 1L)
    val <- val + (alpha / C) * sum(log(pmax(state$class_weights, 1e-300)))
  if (is.null(state$Z_delta))
    for (c in seq_len(C))
      val <- val + (alpha / (C * K)) *
        sum(log(pmax(state$delta_c[[c]], 1e-300)))
  if (Tn > 1L && is.null(state$Z_tau)) {
    n_pat <- K * C
    for (c in seq_len(C)) {
      mats <- state$tau_c[[c]]
      # Under `transition_invariance = "full"` one matrix is estimated and
      # then broadcast, so its prior is carried once, not Tn - 1 times.
      idx  <- if (isTRUE(state$tau_homogeneous)) 1L else seq_along(mats)
      allowed <- state$tau_allowed_c[[c]]
      for (i in idx) for (k in seq_len(K)) {
        ok <- allowed[[i]][k, ]
        if (!any(ok)) next
        val <- val + (alpha / (n_pat * sum(ok))) *
          sum(log(pmax(mats[[i]][k, ok], 1e-300)))
      }
    }
  }
  val
}

.lta_log_prior <- function(state, X, alpha, marginals = NULL) {
  if (!is.null(state$ri))
    return(.lta_ri_log_prior(state, X, alpha, marginals))

  fam <- class(state$mm$models[[1]])[1]
  known <- c("bernoulli", "bernoulli_nan",     # marginal-preserving Beta prior
             "gaussian_unit", "gaussian_unit_nan")  # no measurement prior at all
  if (!fam %in% known) return(NA_real_)

  val <- .lta_log_prior_dt(state, alpha)

  K   <- state$n_statuses
  Tn  <- state$n_times
  a_cat <- .bayes_alpha(state$mm$models[[1]], "categorical")
  if (a_cat > 0 && fam %in% c("bernoulli", "bernoulli_nan")) {
    if (is.null(marginals)) marginals <- .lta_prior_marginals(state, X)
    inv <- .lta_invariant_items(state)
    for (j in seq_len(state$n_items)) {
      grps <- if (j %in% inv) list(seq_len(Tn)) else lapply(seq_len(Tn), identity)
      for (g in grps) {
        p  <- pmin(pmax(state$mm$models[[g[1]]]$parameters$pis[, j], 1e-300),
                   1 - 1e-300)
        mj <- marginals[g[1], j]
        val <- val + length(g) * (a_cat / K) *
          sum(mj * log(p) + (1 - mj) * log1p(-p))
      }
    }
  }
  val
}

# The penalty and its gradient, as a function of the current `state`, returned
# together so the two can never be written from different formulae.
.lta_penalty <- function(state, X, layout, alpha) {
  K <- state$n_statuses
  C <- state$n_classes %||% 1L
  # Read the way m_step.bernoulli() reads it -- off the sub-model, through
  # .bayes_alpha()'s own defaulting -- so the penalty here and the prior in the
  # M-step cannot be two different numbers.
  a_cat <- .bayes_alpha(state$mm$models[[1]], "categorical")
  a_var <- .bayes_alpha(state$mm$models[[1]], "variances")
  val <- 0
  grad <- numeric(0)
  for (b in layout) {
    g <- numeric(b$len)
    if (b$kind == "class" && alpha > 0) {
      # The class weights are themselves a conditional table -- C patterns,
      # one observation apiece -- so they carry the same (alpha / C) form as
      # delta and tau below; see .lta_log_prior(), which this mirrors exactly.
      p <- pmax(state$class_weights, 1e-300)
      val <- val + (alpha / C) * sum(log(p))
      g <- (alpha / C) * (1 - C * p[seq_len(C - 1L)])
    } else if (b$kind == "delta" && alpha > 0) {
      p <- pmax(state$delta_c[[b$c]], 1e-300)
      # Tied: one shared table, so it carries the whole alpha/K rather than a
      # C-th of it -- matching .lta_normalise(patterns = 1L) in the tied M-step,
      # and .lta_log_prior_dt(), which adds the C-th C times over identical
      # vectors and lands on the same number. Untied: the initial-status
      # distribution is estimated once per class, so its (alpha / K) mass is
      # shared C ways, matching .lta_normalise(patterns = C) in the M-step and
      # .lta_log_prior()'s (alpha / (C * K)) term. At C == 1 this is (alpha / K)
      # either way.
      a_d <- if (isTRUE(b$tied)) alpha / K else alpha / (C * K)
      val <- val + a_d * sum(log(p))
      g <- a_d * (1 - K * p[seq_len(K - 1L)])
    } else if (b$kind == "tau" && alpha > 0) {
      p  <- pmax(state$tau_c[[b$c]][[b$i_mat]][b$k, b$allowed], 1e-300)
      m  <- length(p)
      # `alpha` is the whole transition table's mass, spread over its K origin
      # patterns -- times C, when a latent class conditions the table too --
      # matching .lta_normalise(patterns = K * C) in the M-step.
      a_tau <- alpha / (K * C * m)
      val <- val + a_tau * sum(log(p))
      g <- a_tau * (1 - m * p[-m])
    } else if (b$kind == "rho" && a_cat > 0) {
      p  <- pmin(pmax(state$mm$models[[b$grp[1]]]$parameters$pis[, b$j],
                      1e-300), 1 - 1e-300)
      mj <- .lta_rho_prior_marginal(state, X, b)
      # The other program writes one measurement equation per occasion and puts
      # a_cat/K on each, so an item held invariant over the occasions in `grp`
      # carries that mass once per occasion, not once in total. m_step.blocks()
      # applies the same multiplier when it pools the occasions.
      a_rho <- length(b$grp) * a_cat / K
      val <- val + a_rho * sum(mj * log(p) + (1 - mj) * log1p(-p))
      # On the logit scale the derivative collapses to a_rho * (m_j - rho).
      g <- a_rho * (mj - p)
    } else if (b$kind == "log_sd" && a_var > 0) {
      v  <- state$mm$models[[b$grp[1]]]$parameters$covariances[, b$j]
      s2 <- .lta_var_prior_marginal(state, X, b)
      # Truncated inverse-Wishart prior centred on the item's observed
      # marginal variance, spread over `length(grp)` occasions the same way
      # a_rho is above -- matching m_step.gaussian_diag()'s posterior mode
      # sigma^2 = (SS + (alpha/K) s^2) / (n + alpha/K).
      a_v <- length(b$grp) * a_var / K
      val <- val - 0.5 * a_v * sum(log(v) + s2 / v)
      g <- a_v * (s2 / v - 1)
    }
    grad <- c(grad, g)
  }
  list(value = val, gradient = grad)
}

# One L-BFGS climb from wherever EM stopped. Returns `state` unchanged -- never
# a worse fit, and never a fit whose stored posteriors describe other
# parameters -- if the climb does not improve the objective or the model is
# outside the score blocks' scope.
.lta_refine_lbfgs <- function(state, X, alpha = 1.0, max_iter = 200L) {
  if (!.lta_scores_full(state)) return(state)
  w <- state$weights_vec

  layout <- .lta_par_layout(state)
  par0   <- .lta_par_pack(state, layout)

  objective <- function(par) {
    st <- .lta_par_unpack(par, state, layout)
    sc <- .lta_score_matrix(st, X)
    if (is.null(sc)) return(NULL)
    pen <- .lta_penalty(st, X, layout, alpha)
    list(state = st,
         value = sum(w * sc$ll) + pen$value,
         gradient = colSums(sweep(sc$S, 1, w, "*")) + pen$gradient)
  }

  # `optim()` asks for value and gradient separately and asks for both at the
  # same point far more often than not; the forward-backward pass behind them
  # is the expensive part, so it is computed once and reused.
  cache <- NULL
  at <- function(par) {
    if (is.null(cache) || !identical(cache$par, par))
      cache <<- c(list(par = par), objective(par))
    cache
  }

  # Box the logit-scale blocks. Under pure maximum likelihood -- no measurement
  # prior -- the optimum of a sparse LTA can sit ON the boundary, and L-BFGS
  # goes there directly where EM only creeps towards it. That is a real
  # optimum and not a failure, but a rho of exactly 0 or 1 is a degenerate
  # starting point for any later E-step: `refine_from` handed such a donor
  # collapsed by three log-likelihood units, because the responsibilities it
  # implies are themselves degenerate. +/-25 on the logit scale is a
  # probability of 1.4e-11, near enough to the boundary to lose nothing that
  # can be measured and far enough from it to stay a usable fit. Gaussian
  # means are unbounded; they have no boundary to reach.
  bounded <- unlist(lapply(layout, function(b)
    rep(b$kind != "mu", b$len)), use.names = FALSE)
  lo <- ifelse(bounded, -25, -Inf)
  hi <- ifelse(bounded,  25,  Inf)
  # EM can itself arrive at a boundary, which would put the starting vector
  # outside that box; pull it in before optim() is asked to respect it.
  par0 <- pmin(pmax(par0, lo), hi)

  base <- at(par0)
  if (is.null(base$value) || !is.finite(base$value)) return(state)

  opt <- tryCatch(
    stats::optim(par0,
                 fn = function(p) { v <- at(p)$value
                                    if (is.finite(v)) -v else .Machine$double.xmax },
                 gr = function(p) -at(p)$gradient,
                 method = "L-BFGS-B", lower = lo, upper = hi,
                 control = list(maxit = max_iter)),
    error = function(e) NULL)
  if (is.null(opt)) return(state)

  final <- at(opt$par)
  if (!is.finite(final$value) || final$value <= base$value) return(state)

  # The polished parameters need their own E-step before they are returned:
  # gamma, xi and the path entropy all describe the parameters they were
  # computed at, and handing back new parameters with the old posteriors
  # attached is the write-back error d5e1a06 fixed elsewhere in this package.
  # `max_iter = 0` runs no EM iteration at all: .lta_em() falls straight
  # through to its own closing E-step, which is exactly the write-back wanted
  # here -- the polished parameters kept, every posterior recomputed at them.
  out <- .lta_em(final$state, X, max_iter = 0L, alpha = alpha)
  out$refined_lbfgs <- TRUE
  if (out$loglik < state$loglik) return(state)
  # `max_iter = 0` above means .lta_em() ran no iteration, so both of these
  # come back describing the E-step rather than the search: `converged`
  # re-initialised to FALSE and `n_iter` to 0. They belong to the EM run this
  # polish was handed, so they are carried over. Without the second line a
  # polished fit printed "Converged: TRUE (in 0 iterations)".
  out$converged <- state$converged
  out$n_iter    <- state$n_iter
  out
}

# Seed one EM run from a fitted model's own converged parameters. The state is
# built exactly as a random start would build it -- that is what fixes the
# scaffolding the driver expects, the allowed-transition masks and the emission
# skeleton included -- and then every free quantity is overwritten from the
# donor. What is not copied is the regression coefficients: .lta_random_start()
# clears them and the first M-step refits them to the probabilities it is given,
# which are the donor's, so the covariate model resumes at the same point.
.lta_refine_start <- function(state, X, donor) {
  state <- .lta_random_start(state, X)

  state$delta_c <- donor$delta_c
  state$tau_c   <- donor$tau_c
  if (state$n_classes > 1L) state$class_weights <- donor$class_weights

  # An RI fit may take a regular-LTA donor: seed delta, tau and the RI
  # intercepts A from its converged solution and leave the loadings L at their
  # small random draw (roadmap ### 14.10.7). The donor's own `pis` is exactly
  # the K x R intercept table an RI model with L == 0 would report.
  if (!is.null(state$ri) && is.null(donor$ri)) {
    state$ri$A <- qlogis(pmin(pmax(donor$mm$models[[1]]$parameters$pis,
                                   0.05), 0.95))
    return(state)
  }

  mm <- .copy_emission_parameters(state$mm, donor$mm)
  if (is.null(mm))
    stop("`refine_from` has a measurement model of a different shape from the ",
         "one this fit asks for, so its solution cannot be continued.",
         call. = FALSE)
  state$mm <- mm
  state
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

  if (!is.null(state$ri)) {
    J <- state$n_items
    R <- J
    M <- ncol(state$ri$Dnode)
    pis0 <- state$mm$models[[1]]$parameters$pis
    state$ri$A <- qlogis(pmin(pmax(pis0, 0.05), 0.95))
    # Never start a loading at exactly 0: 14.3 says zero is a stationary ridge
    # a restart placed there cannot leave, and the failure it causes -- the fit
    # equals regular LTA -- is indistinguishable from a real "no random
    # intercept needed" result (roadmap ### 14.10.10, failure mode 4).
    state$ri$L <- matrix(stats::runif(R * M, 0.2, 0.8), R, M)
    if (state$ri$kind == "binary")
      state$ri$mass <- rep(1 / length(state$ri$mass), length(state$ri$mass))
    state$mm$models[[1]]$parameters$pis <-
      .lta_ri_integrated_pis(state$ri, K, R)
    for (t in seq_len(Tn))
      state$mm$models[[t]]$parameters$pis <- state$mm$models[[1]]$parameters$pis
  }
  state
}
