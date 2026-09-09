# ==============================================================================
# Direct covariate effects on the indicators (measurement non-invariance)
# ==============================================================================
#
# `predictors_items` asks a different question from every other covariate
# argument fit_lta() takes. `predictors_initial` and `predictors_transition`
# ask which latent status a person is in and where they move; `group` gives
# every group its own prevalences and transitions while the measurement model
# stays invariant. This asks whether the ITEMS mean the same thing to people
# who are in the same latent status -- measurement invariance -- by letting the
# covariate act on each item's linear predictor directly:
#
#   eta[i, k, j, s] = theta[k, j, s] + lambda_j' d_q + x_i' d[k, j, ]
#
# One proportional-odds slope per (latent status, item, covariate), shared
# across occasions. The shift lands in exactly the place the random
# intercept's node shift already lands, which is why the ordinal emission
# needs no new mathematics: see .ordinal_cat_probs() (R/ordinal.R), whose
# K x (S-1) recycling turns a length-K shift into a per-status one.
#
# Two consequences shape everything below.
#
# The covariate enters through a small number of distinct PATTERNS, not one
# value per case. The E-step builds one emission table per (node, pattern) and
# the pattern loop splits the rows rather than multiplying them, so the total
# row-work is unchanged and the cost is loop overhead. That only holds while
# the number of patterns is small, which is why a many-valued covariate is
# refused rather than warned about.
#
# The prior is applied at d = 0. The M-step's grid carries the measurement
# prior on its own rows, with a zero covariate design, so .lta_ri_log_prior()
# and .lta_penalty()'s threshold branches remain correct exactly as written --
# they loop the nodes with no DIF term, which is precisely the rows fitted
# here. Spreading the prior across the covariate patterns instead would force
# four prior functions to change in lockstep, to buy a prior that is set to
# zero in every model this was built for.

# ------------------------------------------------------------------------------
# The design, and the pattern index
# ------------------------------------------------------------------------------

# Built through .lta_ri_design() (R/lta_ri_covariates.R): no intercept column
# and a constant column refused, both of which are right here for the same
# reason they are right there -- a shift shared by every case is absorbed
# exactly by the free per-status thresholds and is not identified.
.lta_dif_init <- function(predictors, n, K, R, label = "predictors_items") {
  Z <- .lta_ri_design(predictors, n, label)

  # match() against unique() rather than factor(): the codes only have to
  # separate the distinct values, and a factor's level ordering is a needless
  # hazard for a numeric covariate. `rows` is built with which() for the same
  # reason -- split()'s ordering is the factor's, not the data's.
  codes <- lapply(seq_len(ncol(Z)), function(j) match(Z[, j], unique(Z[, j])))
  key   <- do.call(paste, c(codes, sep = "\r"))
  u     <- !duplicated(key)
  idx   <- match(key, key[u])
  P     <- sum(u)

  cap <- getOption("mixtureEM.dif_max_patterns", 16L)
  if (P > cap)
    stop(sprintf(paste0(
      "`%s` produced %d distinct covariate patterns, above the limit of %d. ",
      "Each pattern costs one extra emission table per quadrature node in ",
      "every E-step, so a continuous covariate turns a forty-minute fit into ",
      "a multi-day one. Bin or dichotomise the covariate, or raise ",
      "`options(mixtureEM.dif_max_patterns = )` deliberately."),
      label, P, cap), call. = FALSE)

  list(Z = Z, Zu = Z[u, , drop = FALSE], pat = idx, P = P,
       rows = lapply(seq_len(P), function(p) which(idx == p)),
       beta = array(0, c(K, R, ncol(Z))), names = colnames(Z))
}

# The K x R additive shift to every item's linear predictor for covariate
# pattern p. Same scale as the loading's node shift, so it enters at the same
# place: eta[k, j] gains x_p' d[k, j, ].
.lta_dif_shift <- function(dif, p) {
  K <- dim(dif$beta)[1]
  R <- dim(dif$beta)[2]
  s <- matrix(0, K, R)
  for (d in seq_len(ncol(dif$Zu)))
    s <- s + matrix(dif$beta[, , d], K, R) * dif$Zu[p, d]
  s
}

# ------------------------------------------------------------------------------
# The measurement M-step
# ------------------------------------------------------------------------------
#
# A sibling of .lta_ri_mstep_ordinal() (R/lta_ri.R), not a rewrite of it: that
# function is the shipped path and stays untouched, and the duplication here
# buys it zero risk. The difference is one dimension. There, the sufficient
# statistic is K x Q x S_j and the threshold cycle sees one shift per node;
# here it is K x G x S_j over a grid of G = Q*P + Q rows, where a data row is a
# (node, covariate pattern) pair and the last Q rows carry the prior at a zero
# covariate design.
#
# Two cycles, as there. Cycle 1 fits this item's thresholds AND its DIF slopes
# together, one status at a time -- both are status- and item-specific, so
# splitting them would be an extra cycle for nothing. Cycle 2 fits the loading,
# pooling every status and grid row, and is skipped outright when the loading
# is not free (the degenerate one-node factor a no-random-intercept DIF fit
# borrows). Both cycles only increase the penalised objective, so ECM
# monotonicity holds by the same argument as there.
.lta_dif_mstep_ordinal <- function(state, X, E, alpha) {
  K    <- state$n_statuses
  Tn   <- state$n_times
  ri   <- state$ri
  dif  <- state$dif
  Q    <- length(ri$mass)
  P    <- dif$P
  D    <- ncol(dif$Zu)
  w    <- state$weights_vec
  G    <- E$ri$G
  cats <- ri$cats %||% state$mm$models[[1]]$cats
  R    <- length(cats)
  free_loading <- .lta_ri_loading_free(state)

  a_cat     <- .bayes_alpha(state$mm$models[[1]], "categorical")
  prior_obs <- Tn * a_cat / (K * Q)

  # The grid. Data rows first, node-slow and pattern-fast, then the prior's Q
  # rows. `node_of` says which node each row integrates at; `Zg` is the row's
  # covariate design, zero on the prior rows -- that zero is the whole reason
  # the prior functions elsewhere need no DIF term.
  n_grid  <- Q * P + Q
  node_of <- c(rep(seq_len(Q), each = P), seq_len(Q))
  Zg      <- rbind(dif$Zu[rep(seq_len(P), times = Q), , drop = FALSE],
                   matrix(0, Q, D))

  theta  <- ri$theta
  lambda <- ri$L
  beta   <- dif$beta

  for (j in seq_len(R)) {
    Sj   <- cats[j]
    cols <- .ordinal_theta_cols(cats, j)

    # n[k, g, s], accumulated exactly as the RI branch accumulates its
    # K x Q x Sj array but with each occasion loop restricted to the rows of
    # one covariate pattern.
    n_arr <- array(0, dim = c(K, n_grid, Sj))
    for (t in seq_len(Tn)) {
      xj <- X[, .time_block_cols(t, R)[j]]
      for (p in seq_len(P)) {
        rows <- dif$rows[[p]]
        xp   <- xj[rows]
        wp   <- w[rows]
        obs  <- !is.na(xp)
        for (q in seq_len(Q)) {
          g   <- (q - 1L) * P + p            # node-slow, pattern-fast: node_of/Zg
          Gqt <- G[[q]][[t]][rows, , drop = FALSE]
          for (s in seq_len(Sj)) {
            in_s <- obs & (xp == s)
            if (any(in_s))
              n_arr[, g, s] <- n_arr[, g, s] +
                colSums(wp[in_s] * Gqt[in_s, , drop = FALSE])
          }
        }
      }
    }

    # The prior's own rows: the same pseudo-count mass the RI branch uses,
    # spread across this item's categories by the weighted observed marginal,
    # placed where the covariate design is zero.
    m_j <- .lta_ri_item_marginal_ordinal(X, w, R, Tn, j, Sj)
    for (s in seq_len(Sj))
      n_arr[, Q * P + seq_len(Q), s] <-
        n_arr[, Q * P + seq_len(Q), s] + prior_obs * m_j[s]

    theta_j  <- theta[, cols, drop = FALSE]
    lambda_j <- lambda[j, ]
    d_j      <- matrix(beta[, j, ], K, D)
    node_sh  <- as.vector(ri$Dnode %*% lambda_j)[node_of]   # length n_grid

    # Cycle 1: thresholds and DIF slopes together, one status at a time.
    for (k in seq_len(K)) {
      par0 <- c(theta_j[k, ], d_j[k, ])
      nll <- function(par) {
        th_k <- matrix(par[seq_len(Sj - 1L)], 1L, Sj - 1L)
        sh   <- matrix(node_sh + as.vector(Zg %*% par[Sj - 1L + seq_len(D)]),
                       1L, n_grid)
        .ordinal_grid_negloglik(th_k, sh, array(n_arr[k, , ], c(1L, n_grid, Sj)),
                                Sj)
      }
      fit_k <- stats::nlminb(par0, nll)
      theta_j[k, ] <- fit_k$par[seq_len(Sj - 1L)]
      d_j[k, ]     <- fit_k$par[Sj - 1L + seq_len(D)]
    }
    theta[, cols] <- theta_j
    beta[, j, ]   <- d_j

    # Cycle 2: the loading, thresholds and slopes fixed, pooling every status
    # and grid row. Nothing to fit when the factor is the degenerate one.
    if (free_loading) {
      dif_sh <- t(tcrossprod(d_j, Zg))                       # n_grid x K
      nll_l <- function(lam) {
        sh <- t(matrix(as.vector(ri$Dnode %*% lam)[node_of], n_grid, K) + dif_sh)
        .ordinal_grid_negloglik(theta_j, sh, n_arr, Sj)
      }
      lambda[j, ] <- stats::nlminb(lambda_j, nll_l)$par
    }
  }

  ri$theta   <- theta
  ri$L       <- lambda
  ri$cats    <- cats
  dif$beta   <- beta
  state$ri   <- ri
  state$dif  <- dif

  # Node masses: identical to the RI branch, family-agnostic. A DIF fit's
  # factor is either the fixed Gauss-Hermite grid or the degenerate single
  # node, so neither is estimated -- but the branch is kept so the two M-steps
  # stay readable side by side.
  if (identical(ri$kind, "binary")) {
    mass_counts <- vapply(seq_len(Q), function(q) sum(w * E$ri$pq[, q]), numeric(1))
    ri$mass  <- .lta_normalise(mass_counts, alpha)
    state$ri <- ri
  }

  # The reported category probabilities are those of a case with every
  # `predictors_items` covariate at ZERO -- the same convention any regression
  # uses for a "reported" fitted value, and the one lta_covariate_summary()
  # states in its own header. .ordinal_pis_from_theta() is called exactly as
  # the RI branch calls it, with no DIF term.
  state$mm$models[[1]]$parameters$pis <-
    .ordinal_pis_from_theta(theta, lambda, ri$Dnode, ri$mass, cats)
  for (t in seq_len(Tn))
    state$mm$models[[t]]$parameters$pis <- state$mm$models[[1]]$parameters$pis
  state
}
