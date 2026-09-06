# ==============================================================================
# Random-intercept LTA (Muthen & Asparouhov, 2022)
# ==============================================================================
#
# The only thing a random intercept changes about fit_lta() is the emission:
# instead of one K x R table of response logits, each of Q "nodes" (quadrature
# points for a continuous factor, or discrete classes for the binary variant)
# has its own table, built from a shared K x R intercept `A` and an R x M
# loading matrix `L`. See internal/ROADMAP.md ### 14.10 for the full design;
# comments here point back to it rather than re-deriving it.

# state$ri's initial shape, built once in fit_lta() before any EM runs. `A`
# and `L` are NULL until .lta_random_start() draws them.
.lta_ri_init <- function(random_intercept, n_quadrature, n_ri) {
  if (random_intercept == "none") return(NULL)
  if (random_intercept == "continuous") {
    gh <- .gauss_hermite(n_quadrature)
    list(kind = "continuous", Dnode = matrix(gh$z, ncol = 1L), mass = gh$w,
         A = NULL, L = NULL)
  } else {
    Q <- as.integer(n_ri)
    M <- Q - 1L
    list(kind = "binary", Dnode = rbind(0, diag(M)), mass = rep(1 / Q, Q),
         A = NULL, L = NULL)
  }
}

# The numerically-integrated conditional response probability (14.4): what
# gets written into state$mm$models[[t]]$parameters$pis for every downstream
# reader (print(), item_probabilities(), n_parameters(), ...). This has no
# closed form; the published article's eq.-(4) probit approximation is a
# different, less accurate quantity and is never written here.
.lta_ri_integrated_pis <- function(ri, K, R) {
  if (!is.null(ri$theta))
    return(.ordinal_pis_from_theta(ri$theta, ri$L, ri$Dnode, ri$mass, ri$cats))
  Q <- length(ri$mass)
  pis <- matrix(0, K, R)
  for (q in seq_len(Q))
    pis <- pis + ri$mass[q] *
      plogis(ri$A + matrix(ri$L %*% ri$Dnode[q, ], K, R, byrow = TRUE))
  pis
}

# One node's emission table, dispatched on measurement family: ordinal
# (`ri$theta` set) or binary (`ri$A` set). Used by both .lta_ri_e_step() and
# .lta_ri_ll_case(), which otherwise build this identically.
.lta_ri_node_pis <- function(ri, K, R, q) {
  if (!is.null(ri$theta))
    return(.ordinal_node_pis(ri$theta, ri$L, ri$Dnode, ri$cats, q))
  plogis(ri$A + matrix(ri$L %*% ri$Dnode[q, ], K, R, byrow = TRUE))
}

# The Q emission tables shared by every RI likelihood computation: one node's
# response-probability table, built at every case. Factored out of
# .lta_ri_e_step() and .lta_ri_ll_case() so the two can never build it
# differently -- they are asserted to agree to 1e-10
# (tests/testthat/test-lta-ri.R:577) and the finite-difference Hessian
# differentiates only the second.
.lta_ri_logB <- function(state, X) {
  K  <- state$n_statuses
  Tn <- state$n_times
  R  <- state$n_items
  ri <- state$ri
  Q  <- length(ri$mass)
  logB <- vector("list", Q)
  for (q in seq_len(Q)) {
    pis_q <- .lta_ri_node_pis(ri, K, R, q)
    mm_q  <- state$mm
    for (t in seq_len(Tn)) mm_q$models[[t]]$parameters$pis <- pis_q
    logB[[q]] <- .lta_emission_loglik(mm_q, X)
  }
  logB
}

# Weighted observed marginal of item j, pooled over every occasion -- RI fits
# require measurement_invariance = "full" (fit_lta() refuses otherwise), so
# "pooled over the occasions this item's M-step pools" is always all of them.
# Mirrors .lta_rho_prior_marginal()'s treatment of missing data exactly.
.lta_ri_item_marginal <- function(X, w, R, Tn, j) {
  num <- 0; den <- 0
  for (t in seq_len(Tn)) {
    xj  <- X[, .time_block_cols(t, R)[j]]
    obs <- !is.na(xj)
    num <- num + sum(w[obs] * xj[obs])
    den <- den + sum(w[obs])
  }
  if (den == 0) 0.5 else num / den
}

# ------------------------------------------------------------------------------
# E-step: Q forward-backward passes per class, mixed by the node posterior,
# then (with more than one class) by the class posterior
# ------------------------------------------------------------------------------
#
# Returns the shape .lta_em()'s own e_step() returns, so nothing downstream
# (the per-class delta/tau M-steps, .lta_pair_counts(), .lta_mixed_gamma(),
# state$abs_ent_path) needs to know a random intercept is there at all. The
# one extra field, `ri = list(pq, G)`, is read only by .lta_ri_mstep() and by
# the sign-normalisation step, and carries no class dimension: the factor
# (`A`, `L`, and the binary variant's `mass`) is shared across classes -- a
# person's response tendency doesn't depend on which chain they follow -- so
# its sufficient statistic is the (status, node) posterior pooled over class,
# weighted by each case's own class posterior (14.7 Phase B: "a double outer
# loop over class and node, and no new mathematics").
.lta_ri_e_step <- function(state, X, w) {
  n  <- nrow(X)
  K  <- state$n_statuses
  Tn <- state$n_times
  R  <- state$n_items
  ri <- state$ri
  Q  <- length(ri$mass)
  C  <- state$n_classes %||% 1L

  # The emission table at each node depends only on the shared factor, not on
  # class, so it is built once per node and reused across every class's
  # forward-backward pass below.
  logB  <- .lta_ri_logB(state, X)
  lmass <- .lta_ri_log_mass(state, n)

  # Inner mix, per class: exactly the single-class node loop this function
  # has always run, using that class's own delta/tau.
  fb   <- vector("list", C)
  ll_c <- matrix(0, n, C)
  pq_c <- vector("list", C)
  for (c in seq_len(C)) {
    sub <- .lta_class_state(state, c)
    log_delta <- .lta_log_delta(sub)
    log_tau   <- .lta_log_tau(sub)

    # keep_pairwise = TRUE always: the pairwise blocks must be re-weighted by
    # the node posterior before the transition M-step sees them (14.10.10,
    # failure mode 5 if this is ever made conditional).
    fb[[c]] <- lapply(seq_len(Q), function(q)
      .lta_forward_backward(logB[[q]], log_delta, log_tau, w,
                            keep_pairwise = TRUE))

    lq <- vapply(seq_len(Q), function(q)
      lmass[, q] + fb[[c]][[q]]$ll, numeric(n))
    ll_c[, c] <- logsumexp(lq, MARGIN = 1)
    pq_c[[c]] <- exp(lq - ll_c[, c])
  }

  # Outer mix, over class -- identical in shape to the plain (no-RI) e_step's
  # own C > 1 branch in .lta_em(), just fed `ll_c` in place of a directly
  # computed per-class `fb$ll`.
  if (C == 1L) {
    ll   <- ll_c[, 1]
    post <- matrix(1, n, 1L)
  } else {
    lp   <- sweep(ll_c, 2, log(pmax(state$class_weights, 1e-300)), "+")
    ll   <- logsumexp(lp, MARGIN = 1)
    post <- exp(lp - ll)
  }

  # Per-class, node-mixed status/pairwise posteriors, and the per-class JOINT
  # (status, node) posterior computed before it is ever summed over q. This
  # is the single easiest thing in the Part to get wrong (14.10.10, failure
  # mode 1): the measurement M-step needs the joint quantities below, never
  # pq[, q] times an already node-mixed gamma/pairwise.
  es     <- vector("list", C)
  Gjoint <- vector("list", C)
  for (c in seq_len(C)) {
    Gc <- lapply(seq_len(Q), function(q)
      lapply(fb[[c]][[q]]$gamma, function(g) g * pq_c[[c]][, q]))
    Pc <- lapply(seq_len(Q), function(q)
      lapply(fb[[c]][[q]]$pairwise, function(pt)
        lapply(pt, function(block) block * pq_c[[c]][, q])))

    gamma_mix <- lapply(seq_len(Tn), function(t)
      Reduce(`+`, lapply(seq_len(Q), function(q) Gc[[q]][[t]])))
    pair_mix <- if (Tn > 1L) lapply(seq_len(Tn - 1L), function(t)
      lapply(seq_len(K), function(k)
        Reduce(`+`, lapply(seq_len(Q), function(q) Pc[[q]][[t]][[k]])))) else
      list()

    es[[c]]     <- list(ll = ll_c[, c], gamma = gamma_mix, xi = NULL,
                        pairwise = pair_mix)
    Gjoint[[c]] <- Gc
  }

  # The shared factor's sufficient statistic: the per-class joint (status,
  # node) posterior, pooled over class weighted by each case's class
  # posterior. With one class this reduces to Gjoint[[1]] unchanged.
  G <- lapply(seq_len(Q), function(q)
    lapply(seq_len(Tn), function(t) {
      terms <- lapply(seq_len(C), function(c) {
        g <- Gjoint[[c]][[q]][[t]]
        if (C == 1L) g else g * post[, c]
      })
      Reduce(`+`, terms)
    }))
  pq <- if (C == 1L) pq_c[[1]] else
    Reduce(`+`, lapply(seq_len(C), function(c) pq_c[[c]] * post[, c]))

  list(es = es, post = post, ll = ll, ri = list(pq = pq, G = G))
}

# The per-case observed-data log-likelihood of an RI fit: one forward-backward
# per node, mixed over nodes AFTER the chain, then over classes. This is what
# the finite-difference Hessian differentiates, so it must be the same
# quantity .lta_ri_e_step() returns as `ll` -- the RI packing test asserts
# that to 1e-10.
#
# `.lta_ri_e_step()` itself is deliberately not reused: it runs with
# `keep_pairwise = TRUE` and builds the joint (status, node) posteriors, and
# the Hessian below calls this 2p(p+1) times and wants none of it. That is not
# an optimisation, it is what makes the sandwich affordable at Q = 50.
.lta_ri_ll_case <- function(state, X) {
  n  <- nrow(X)
  K  <- state$n_statuses
  Tn <- state$n_times
  R  <- state$n_items
  ri <- state$ri
  Q  <- length(ri$mass)
  C  <- state$n_classes %||% 1L
  w  <- state$weights_vec

  logB  <- .lta_ri_logB(state, X)
  lmass <- .lta_ri_log_mass(state, n)

  ll_c <- matrix(0, n, C)
  for (c in seq_len(C)) {
    sub <- .lta_class_state(state, c)
    log_delta <- .lta_log_delta(sub)
    log_tau   <- .lta_log_tau(sub)
    lq <- vapply(seq_len(Q), function(q)
      lmass[, q] +
        .lta_forward_backward(logB[[q]], log_delta, log_tau, w)$ll, numeric(n))
    ll_c[, c] <- logsumexp(lq, MARGIN = 1)
  }
  if (C == 1L) return(ll_c[, 1])
  logsumexp(sweep(ll_c, 2, log(pmax(state$class_weights, 1e-300)), "+"),
            MARGIN = 1)
}

# ------------------------------------------------------------------------------
# Measurement M-step: one aggregated binomial GLM per item
# ------------------------------------------------------------------------------
#
# Replaces the plain m_step(state$mm, ...) call for RI fits. The linear
# predictor depends on the data only through (status k, node q), so pooling
# every case and every occasion into a (K * Q) x 2 table per item is exact,
# not an approximation -- 14.10.5's whole point, and what keeps an RI fit
# affordable at Q = 20.
.lta_ri_mstep <- function(state, X, E, alpha) {
  if (!is.null(state$ri$theta) || .lta_is_ordinal_model(state$mm$models[[1]]))
    return(.lta_ri_mstep_ordinal(state, X, E, alpha))

  K  <- state$n_statuses
  R  <- state$n_items
  Tn <- state$n_times
  ri <- state$ri
  Q  <- length(ri$mass)
  M  <- ncol(ri$Dnode)
  w  <- state$weights_vec
  G  <- E$ri$G

  # Design: K status dummies (k fastest) then M node columns, expand.grid
  # order -- no separate intercept column, the K dummies already saturate it.
  D <- cbind(diag(K)[rep(seq_len(K), times = Q), , drop = FALSE],
             ri$Dnode[rep(seq_len(Q), each = K), , drop = FALSE])

  a_cat <- .bayes_alpha(state$mm$models[[1]], "categorical")
  prior_obs <- Tn * a_cat / (K * Q)

  A <- ri$A
  L <- ri$L
  for (j in seq_len(R)) {
    succ <- matrix(0, K, Q)
    tot  <- matrix(0, K, Q)
    for (t in seq_len(Tn)) {
      xj  <- X[, .time_block_cols(t, R)[j]]
      obs <- !is.na(xj)
      for (q in seq_len(Q)) {
        Gqt <- G[[q]][[t]]
        succ[, q] <- succ[, q] + colSums(w[obs] * Gqt[obs, , drop = FALSE] * xj[obs])
        tot[, q]  <- tot[, q]  + colSums(w[obs] * Gqt[obs, , drop = FALSE])
      }
    }
    # as.vector() on a K x Q matrix is column-major -- k fastest, q slowest --
    # which is exactly D's row order above.
    succ_vec <- as.vector(succ)
    tot_vec  <- as.vector(tot)

    m_j <- .lta_ri_item_marginal(X, w, R, Tn, j)
    pri <- .wglm_prior_rows(D, rep(m_j, K * Q), prior_obs)

    fit_j <- .wglm_fit(rbind(D, pri$D),
                       y = c(succ_vec / tot_vec, pri$y),
                       w = c(tot_vec, pri$w),
                       fam = .wglm_family("binomial"),
                       start = c(A[, j], L[j, ]))
    A[, j] <- fit_j$coefficients[seq_len(K)]
    L[j, ] <- fit_j$coefficients[K + seq_len(M)]
  }
  ri$A <- A
  ri$L <- L

  # The binary variant's node masses are themselves estimated; the continuous
  # variant's Gauss-Hermite weights are fixed and must never be touched here.
  if (ri$kind == "binary") {
    mass_counts <- vapply(seq_len(Q), function(q) sum(w * E$ri$pq[, q]), numeric(1))
    ri$mass <- .lta_normalise(mass_counts, alpha)
  }

  state$ri <- ri
  state$mm$models[[1]]$parameters$pis <- .lta_ri_integrated_pis(ri, K, R)
  for (t in seq_len(Tn))
    state$mm$models[[t]]$parameters$pis <- state$mm$models[[1]]$parameters$pis
  state
}

# ------------------------------------------------------------------------------
# Measurement M-step, ordinal family: two-cycle ECM, one item at a time
# ------------------------------------------------------------------------------
#
# The binomial GLM above collapses each item to one aggregated (K*Q) x 2
# table because a Bernoulli response has one sufficient statistic. An ordinal
# item's cumulative-logit response does not collapse the same way -- its
# sufficient statistic is a K x Q x S_r table of category counts (### 14.18,
# W4) -- so this is not `.wglm_fit()` on a wider table, it is a small
# per-item Newton problem: thresholds first (lambda fixed, one status at a
# time), then the loading (thresholds fixed, pooling every status and node).
# Both cycles only ever increase the penalised objective, so alternating them
# preserves EM monotonicity exactly as alternating the delta/tau cycles does
# elsewhere in this file.
.lta_ri_mstep_ordinal <- function(state, X, E, alpha) {
  K    <- state$n_statuses
  Tn   <- state$n_times
  ri   <- state$ri
  Q    <- length(ri$mass)
  w    <- state$weights_vec
  G    <- E$ri$G
  cats <- ri$cats %||% state$mm$models[[1]]$cats
  R    <- length(cats)

  a_cat     <- .bayes_alpha(state$mm$models[[1]], "categorical")
  prior_obs <- Tn * a_cat / (K * Q)

  theta  <- ri$theta
  lambda <- ri$L

  for (j in seq_len(R)) {
    Sj   <- cats[j]
    cols <- .ordinal_theta_cols(cats, j)

    # n[k, q, s]: the K x Q x Sj sufficient statistic, accumulated from the
    # joint (status, node) posterior E$ri$G exactly as the binary branch
    # accumulates succ/tot, but keeping every category separate.
    n_arr <- array(0, dim = c(K, Q, Sj))
    for (t in seq_len(Tn)) {
      xj  <- X[, .time_block_cols(t, R)[j]]
      obs <- !is.na(xj)
      for (q in seq_len(Q)) {
        Gqt <- G[[q]][[t]]
        for (s in seq_len(Sj)) {
          in_s <- obs & (xj == s)
          if (any(in_s))
            n_arr[, q, s] <- n_arr[, q, s] +
              colSums(w[in_s] * Gqt[in_s, , drop = FALSE])
        }
      }
    }

    # Dirichlet-style pseudo-count prior: the same prior_obs mass the binary
    # branch uses, spread across this item's S_r categories by the weighted
    # observed marginal (mirrors m_step.ordinal's marginal_prob, R/ordinal.R).
    m_j <- .lta_ri_item_marginal_ordinal(X, w, R, Tn, j, Sj)
    for (s in seq_len(Sj))
      n_arr[, , s] <- n_arr[, , s] + prior_obs * m_j[s]

    theta_j  <- theta[, cols, drop = FALSE]
    lambda_j <- lambda[j, ]

    # Cycle 1: per-status Newton on thresholds, lambda fixed. matrix() below
    # guards against R silently dropping the Q dimension when Q == 1 (the
    # n_quadrature = 1 test case).
    shift <- as.vector(ri$Dnode %*% lambda_j)
    for (k in seq_len(K))
      theta_j[k, ] <- .ordinal_newton_theta(
        theta_j[k, ], matrix(n_arr[k, , ], Q, Sj), shift, Sj)
    theta[, cols] <- theta_j

    # Cycle 2: one Newton step on the loading, thresholds fixed, pooling
    # every status and node.
    lambda[j, ] <- .ordinal_newton_lambda(lambda_j, theta_j, ri$Dnode, n_arr, Sj)
  }
  ri$theta <- theta
  ri$L     <- lambda
  ri$cats  <- cats

  # Node masses: identical to the binary branch, family-agnostic.
  if (ri$kind == "binary") {
    mass_counts <- vapply(seq_len(Q), function(q) sum(w * E$ri$pq[, q]), numeric(1))
    ri$mass <- .lta_normalise(mass_counts, alpha)
  }

  state$ri <- ri
  state$mm$models[[1]]$parameters$pis <-
    .ordinal_pis_from_theta(theta, lambda, ri$Dnode, ri$mass, cats)
  for (t in seq_len(Tn))
    state$mm$models[[t]]$parameters$pis <- state$mm$models[[1]]$parameters$pis
  state
}

# The measurement side of .lta_log_prior() for an RI fit: the closed form of
# the same pseudo-observation prior .lta_ri_mstep() fits against, so the EM
# stopping rule and the restart ranking climb the objective the M-step
# actually climbs (see ### 41.1 re-opened in internal/ROADMAP.md for what goes
# wrong when they do not match). The delta/tau terms are unaffected by the
# measurement family and are shared with the regular-LTA path via
# .lta_log_prior_dt().
.lta_ri_log_prior <- function(state, X, alpha, marginals = NULL) {
  val <- .lta_log_prior_dt(state, alpha)

  ri <- state$ri
  K  <- state$n_statuses
  R  <- state$n_items
  Tn <- state$n_times
  Q  <- length(ri$mass)
  a_cat <- .bayes_alpha(state$mm$models[[1]], "categorical")
  if (a_cat > 0 && !is.null(ri$theta)) {
    w <- state$weights_vec
    prior_obs <- Tn * a_cat / (K * Q)
    cats <- ri$cats
    for (j in seq_len(R)) {
      Sj    <- cats[j]
      m_j   <- .lta_ri_item_marginal_ordinal(X, w, R, Tn, j, Sj)
      cols  <- .ordinal_theta_cols(cats, j)
      theta_j <- ri$theta[, cols, drop = FALSE]
      for (q in seq_len(Q)) {
        shift <- sum(ri$L[j, ] * ri$Dnode[q, ])
        p <- pmin(pmax(.ordinal_cat_probs(theta_j, shift, Sj), 1e-300), 1 - 1e-300)
        val <- val + prior_obs * sum(sweep(log(p), 2, m_j, "*"))
      }
    }
  } else if (a_cat > 0) {
    w <- state$weights_vec
    prior_obs <- Tn * a_cat / (K * Q)
    for (j in seq_len(R)) {
      m_j <- .lta_ri_item_marginal(X, w, R, Tn, j)
      for (q in seq_len(Q)) {
        eta <- ri$A[, j] + drop(ri$Dnode[q, , drop = FALSE] %*% ri$L[j, ])
        p <- pmin(pmax(plogis(eta), 1e-300), 1 - 1e-300)
        val <- val + prior_obs * sum(m_j * log(p) + (1 - m_j) * log1p(-p))
      }
    }
  }
  # The binary variant's node masses are themselves estimated (the
  # continuous variant's Gauss-Hermite weights are fixed and carry no prior
  # at all), by .lta_normalise(mass_counts, alpha) at .lta_ri_mstep() above --
  # an implicit Dirichlet(alpha/Q, ..., alpha/Q), the same (alpha / patterns)
  # form as the class-mixing term in .lta_log_prior_dt() with Q in place of C.
  # Missing this term here would let the M-step's actual prior and the
  # quantity EM's monotonicity check (and the L-BFGS polish, .lta_penalty())
  # climb disagree on the one place they must not.
  if (alpha > 0 && identical(ri$kind, "binary") && Q > 1L)
    val <- val + (alpha / Q) * sum(log(pmax(ri$mass, 1e-300)))
  val
}

# ------------------------------------------------------------------------------
# Sign / label normalisation, run once at the end of every .lta_em() call
# ------------------------------------------------------------------------------
#
# A random intercept is identified only up to a reflection (continuous) or a
# relabelling of which node is the anchor (binary) -- the likelihood is
# exactly unchanged by either. Run after every EM call, not only the winning
# restart: the reparametrisation preserves score_of() exactly, so normalising
# early costs nothing and means a restart's stored `ri$pq` (once recomputed by
# the following e_step) never needs its own separate permutation.
#
# Continuous: flip on the loading with the largest absolute value, not the
# first -- on a fixture where several loadings are indistinguishable from
# zero, normalising on a near-zero one makes the reported sign a coin flip
# between restarts (14.10.7). Binary: relabel so the anchored (zero) node is
# the one with the largest mass, carrying `mass` and the columns of `L` with
# it; A absorbs the old anchor's shift so every node's actual per-item logit
# (A[k,j] + Dnode[q,] %*% L[j,]) is unchanged, node for node.
.lta_ri_sign_normalise <- function(state) {
  ri <- state$ri
  K  <- state$n_statuses
  R  <- state$n_items
  A  <- ri$A

  if (ri$kind == "continuous") {
    j_star <- which.max(abs(ri$L[, 1]))
    if (ri$L[j_star, 1] < 0) {
      ri$L     <- -ri$L
      ri$Dnode <- -ri$Dnode   # symmetric nodes: the integral is unchanged
      # promote() (R/lta.R) later overwrites Dnode with the unflipped full
      # grid, licensed only by L, Dnode and ri_beta moving together -- flip
      # two of the three and it silently changes the model.
      if (!is.null(state$ri_beta)) state$ri_beta <- -state$ri_beta
    }
  } else {
    Q <- length(ri$mass)
    delta <- matrix(0, Q, R)
    for (q in seq_len(Q))
      delta[q, ] <- as.vector(ri$Dnode[q, , drop = FALSE] %*% t(ri$L))
    q_star <- which.max(ri$mass)
    if (q_star != 1L) {
      perm      <- c(q_star, setdiff(seq_len(Q), q_star))
      new_delta <- sweep(delta[perm, , drop = FALSE], 2, delta[q_star, ], "-")
      if (!is.null(ri$theta)) {
        # Shifting item j's whole linear predictor by a constant is the same
        # as shifting only the base threshold theta_2: the increments
        # (eta_3, ..., eta_S) are relative gaps and are untouched.
        for (j in seq_len(R)) {
          base_col <- .ordinal_theta_cols(ri$cats, j)[1]
          ri$theta[, base_col] <- ri$theta[, base_col] + delta[q_star, j]
        }
      } else {
        A <- A + matrix(delta[q_star, ], K, R, byrow = TRUE)
      }
      ri$L    <- t(new_delta[-1, , drop = FALSE])
      ri$mass <- ri$mass[perm]
    }
  }

  if (any(abs(ri$L) > 10))
    warning("A random-intercept loading exceeds 10 in absolute value; the ",
            "model may be weakly identified on these data.", call. = FALSE)

  if (is.null(ri$theta)) ri$A <- A
  state$ri <- ri
  state$mm$models[[1]]$parameters$pis <- .lta_ri_integrated_pis(ri, K, R)
  for (t in seq_len(state$n_times))
    state$mm$models[[t]]$parameters$pis <- state$mm$models[[1]]$parameters$pis
  state
}

#' Random-intercept factor scores
#'
#' The between-subject factor a random-intercept `fit_lta()` model estimates,
#' one row per case. `map_class` is the node at the posterior mode, and every
#' node's own posterior probability is returned as `node1_posterior`,
#' `node2_posterior`, and so on. For `random_intercept = "continuous"`, where
#' the nodes are quadrature points on a single scale rather than unordered
#' classes, `mean_score` (the posterior mean) and `map_score` (the quadrature
#' node at the posterior mode) are added as well. When the fit used
#' `predictors_random_intercept`, `predicted_mean` (the regression-implied
#' factor mean, `x_i'beta`) is added too, alongside the posterior `mean_score`.
#'
#' @param fit A model fitted by [`fit_lta()`] with `random_intercept` not
#'   `"none"`.
#' @return A data frame with one row per case in the fitted data.
#' @export
random_intercept_scores <- function(fit) {
  if (is.null(fit$ri) || is.null(fit$ri$pq))
    stop("`fit` was not fitted with `random_intercept`.", call. = FALSE)
  pq <- fit$ri$pq
  z  <- if (ncol(fit$ri$Dnode) == 1L) fit$ri$Dnode[, 1] else NULL
  map_idx <- max.col(pq, ties.method = "first")
  out <- data.frame(map_class = map_idx)
  if (!is.null(z)) {
    out$mean_score <- as.vector(pq %*% z)
    out$map_score  <- z[map_idx]
  }
  if (!is.null(fit$Z_ri) && !is.null(fit$ri_beta))
    out$predicted_mean <- as.vector(fit$Z_ri %*% fit$ri_beta)
  colnames(pq) <- paste0("node", seq_len(ncol(pq)), "_posterior")
  cbind(out, pq)
}

# ------------------------------------------------------------------------------
# Warm start: place the statuses before the factor is switched on
# ------------------------------------------------------------------------------
#
# A random intercept and the latent statuses compete to explain the same
# associations among the items, and from a random start the factor can win --
# EM then settles on a solution with the right statuses but inflated loadings,
# a genuine local maximum a short way below the global one. Measured on the
# LTA-FAQ benchmark (roadmap ### 14.10), fifty random restarts converged to
# -14443.775 against the two reference programs' -14442.017, with loadings of
# 1.0 to 4.5 against their 0.12 to 0.97. Those programs pay for that surface
# with 250 and 1000 restarts respectively.
#
# Paying instead with a better start is what Muthen & Asparouhov (2022)
# themselves recommend, and it is what `refine_from` already does by hand: run
# the model WITHOUT the random intercept first, then turn the factor on from
# there. This does it for each restart automatically, so the search that a user
# gets by default is a search over where the factor starts rather than over the
# statuses and the factor at once. The plain pass is cheap -- one iteration of
# it costs a `Q`th of an RI iteration, since the RI E-step runs one whole
# forward-backward per node.
#
# It is deliberately NOT run to convergence. Every restart would then reach the
# same plain-LTA optimum and the pool would differ only in its loading draw,
# throwing away the diversity the restarts are there to provide.
#
# fit_lta() calls this for the CONTINUOUS variant only, and that restriction is
# measured rather than cautious. The binary variant's nodes are estimated
# classes -- `n_ri` IS part of the model, the same distinction the quadrature
# ladder draws -- so its search is an ordinary mixture search that wants
# diverse starts, and seeding every restart from one plain-LTA solution takes
# that diversity away. On the LTA-FAQ binary benchmark fifty restarts reach the
# reference programs' -14435.618 exactly without this, and stop 3.78 short with
# it; 200 restarts recover it, at four times the cost. The continuous variant
# is the opposite case: this is what closes its 1.76.
.lta_ri_warm_start <- function(state, X, alpha, max_iter = 200L) {
  ri <- state$ri
  plain <- state
  plain$ri      <- NULL
  plain$ri_beta <- NULL
  warm <- try(.lta_em(plain, X, max_iter = max_iter, tol = 1e-6, alpha = alpha),
              silent = TRUE)
  if (inherits(warm, "try-error")) return(state)

  # The plain pass's own `pis` is exactly the table a random intercept with
  # zero loadings would report, which is where the intercept parameters
  # belong; the loadings stay at the draw this restart was given. Mirrors
  # .lta_refine_start()'s treatment of a regular-LTA donor.
  if (!is.null(ri$theta)) {
    cats     <- ri$cats %||% state$mm$models[[1]]$cats
    ri$theta <- .ordinal_theta_from_pis(warm$mm$models[[1]]$parameters$pis, cats)
    ri$cats  <- cats
  } else {
    ri$A <- qlogis(pmin(pmax(warm$mm$models[[1]]$parameters$pis, 0.05), 0.95))
  }
  warm$ri      <- ri
  warm$ri_beta <- state$ri_beta
  warm$mm$models[[1]]$parameters$pis <-
    .lta_ri_integrated_pis(ri, state$n_statuses, state$n_items)
  for (t in seq_len(state$n_times))
    warm$mm$models[[t]]$parameters$pis <-
      warm$mm$models[[1]]$parameters$pis
  warm
}
