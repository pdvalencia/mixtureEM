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
#
# `mode` says which slopes are free: an R x D integer matrix, 0 for none, 1 for
# one slope shared by every status (uniform DIF) and 2 for one per status
# (non-uniform). The matrix form of `predictors_items` is every item by status,
# which is what `mode = NULL` means. `beta` keeps its K x R x D shape whatever
# the mode, so every reader of it stays as written: a uniform slope is K equal
# rows and an absent one is zero.
.lta_dif_init <- function(predictors, n, K, R, label = "predictors_items",
                          mode = NULL) {
  Z <- .lta_ri_design(predictors, n, label)
  if (is.null(mode)) {
    mode <- matrix(2L, R, ncol(Z))
  } else {
    mode <- mode[, colnames(Z), drop = FALSE]
  }
  dimnames(mode) <- NULL

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
       beta = array(0, c(K, R, ncol(Z))), names = colnames(Z), mode = mode)
}

# The list form of `predictors_items`: item name -> the covariates acting on
# that item, as in fit_mixture() (.dif_spec(), R/categorical_dif.R). Items are
# named as the fit names them, which is occasion 1's column names with the
# occasion marker dropped; the slope is shared across occasions either way.
# Returns the union design (one column per distinct covariate, by name) and the
# R x D mode matrix .lta_dif_init() reads. A list entry is uniform unless its
# item is named in `by_status`.
.lta_dif_spec <- function(spec, by_status, item_names, n) {
  if (is.null(names(spec)) || any(!nzchar(names(spec))))
    stop("`predictors_items` given as a list must be named: item name = ",
         "the covariates acting on that item, e.g. list(item3 = ",
         "data.frame(female)).", call. = FALSE)
  bad <- setdiff(c(names(spec), by_status), item_names)
  if (length(bad))
    stop(sprintf("`predictors_items` names %s not among the items (%s): %s.",
                 if (length(bad) > 1L) "items" else "an item",
                 paste(item_names, collapse = ", "),
                 paste(bad, collapse = ", ")), call. = FALSE)
  bad <- setdiff(by_status, names(spec))
  if (length(bad))
    stop(sprintf(paste0("`predictors_items_by_status` names %s with no entry ",
                        "in `predictors_items`: %s."),
                 if (length(bad) > 1L) "items" else "an item",
                 paste(bad, collapse = ", ")), call. = FALSE)
  # Columns are matched across items by name, so a nameless vector could not
  # say whether two items share a covariate or have two of their own.
  unnamed <- vapply(spec, function(x) is.null(colnames(x)), logical(1))
  if (any(unnamed))
    stop(sprintf(paste0("Each entry of `predictors_items` must be a data frame ",
                        "or a matrix with column names (%s has none), e.g. ",
                        "list(item3 = data.frame(female))."),
                 paste(names(spec)[unnamed], collapse = ", ")), call. = FALSE)
  per_item <- lapply(names(spec), function(it)
    as.matrix(prepare_covariates(spec[[it]])))
  if (any(vapply(per_item, nrow, 1L) != n))
    stop("Every entry of `predictors_items` must have one row per case.",
         call. = FALSE)
  covs <- unique(unlist(lapply(per_item, colnames)))
  Z    <- matrix(NA_real_, n, length(covs), dimnames = list(NULL, covs))
  mode <- matrix(0L, length(item_names), length(covs),
                 dimnames = list(item_names, covs))
  for (i in seq_along(per_item)) {
    cn <- colnames(per_item[[i]])
    for (cc in cn) {
      if (!all(is.na(Z[, cc])) &&
          !isTRUE(all.equal(Z[, cc], per_item[[i]][, cc], check.attributes = FALSE)))
        stop(sprintf(paste0("`predictors_items` gives two different columns ",
                            "named `%s`; a covariate must be the same variable ",
                            "on every item it acts on."), cc), call. = FALSE)
      Z[, cc] <- per_item[[i]][, cc]
    }
    mode[names(spec)[i], cn] <- if (names(spec)[i] %in% by_status) 2L else 1L
  }
  list(Z = Z, mode = mode)
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

    # Cycle 1: thresholds and the status-specific DIF slopes together, one
    # status at a time; a uniform or absent slope is held where it is.
    by_st <- which(.lta_dif_mode(dif)[j, ] == 2L)
    unif  <- which(.lta_dif_mode(dif)[j, ] == 1L)
    for (k in seq_len(K)) {
      par0 <- c(theta_j[k, ], d_j[k, by_st])
      nll <- function(par) {
        th_k <- matrix(par[seq_len(Sj - 1L)], 1L, Sj - 1L)
        d_k  <- d_j[k, ]
        d_k[by_st] <- par[Sj - 1L + seq_along(by_st)]
        sh   <- matrix(node_sh + as.vector(Zg %*% d_k), 1L, n_grid)
        .ordinal_grid_negloglik(th_k, sh, array(n_arr[k, , ], c(1L, n_grid, Sj)),
                                Sj)
      }
      fit_k <- stats::nlminb(par0, nll)
      theta_j[k, ]     <- fit_k$par[seq_len(Sj - 1L)]
      d_j[k, by_st]    <- fit_k$par[Sj - 1L + seq_along(by_st)]
    }

    # Cycle 1b: the uniform slopes, one value shared by every status, pooling
    # the statuses with the thresholds held -- the loading cycle's shape. It
    # only increases the objective, so monotonicity is kept.
    if (length(unif)) {
      nll_u <- function(u) {
        d_all <- d_j
        d_all[, unif] <- matrix(u, K, length(unif), byrow = TRUE)
        sh <- t(matrix(node_sh, n_grid, K) + tcrossprod(Zg, d_all))
        .ordinal_grid_negloglik(theta_j, sh, n_arr, Sj)
      }
      u <- stats::nlminb(d_j[1L, unif], nll_u)$par
      d_j[, unif] <- matrix(u, K, length(unif), byrow = TRUE)
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

# ------------------------------------------------------------------------------
# The free slopes of one item, as the parameter vector carries them
# ------------------------------------------------------------------------------
#
# Covariate slowest, status fastest -- as.vector() on the K x D slope matrix
# when every slope is by status, which is the order the matrix form has always
# been packed in. A uniform slope contributes one entry and an absent one none.
# .lta_par_layout(), .lta_par_pack(), .lta_par_unpack() and
# .lta_score_matrix() all go through these, so the four cannot drift apart.
# A fit saved before `mode` existed had every slope by status.
.lta_dif_mode <- function(dif)
  dif$mode %||% matrix(2L, dim(dif$beta)[2], dim(dif$beta)[3])

.lta_dif_len <- function(dif, j) {
  K    <- dim(dif$beta)[1]
  mode <- .lta_dif_mode(dif)
  as.integer(sum(mode[j, ] == 1L) + K * sum(mode[j, ] == 2L))
}

.lta_dif_pack <- function(dif, j) {
  mode <- .lta_dif_mode(dif)
  unlist(lapply(seq_len(ncol(dif$Zu)), function(d)
    switch(mode[j, d] + 1L, NULL, dif$beta[1L, j, d], dif$beta[, j, d])),
    use.names = FALSE)
}

.lta_dif_unpack <- function(dif, j, v) {
  K   <- dim(dif$beta)[1]
  pos <- 0L
  for (d in seq_len(ncol(dif$Zu))) {
    m <- .lta_dif_mode(dif)[j, d]
    if (m == 1L) {
      dif$beta[, j, d] <- v[pos + 1L]; pos <- pos + 1L
    } else if (m == 2L) {
      dif$beta[, j, d] <- v[pos + seq_len(K)]; pos <- pos + K
    }
  }
  dif
}

# `S` holds one column per (covariate, status), status fastest, which is what a
# by-status slope is scored on. A uniform slope moves every status's linear
# predictor at once, so its score is the sum of its K columns.
.lta_dif_score <- function(dif, j, S) {
  K <- dim(dif$beta)[1]
  do.call(cbind, lapply(seq_len(ncol(dif$Zu)), function(d) {
    cols <- (d - 1L) * K + seq_len(K)
    switch(.lta_dif_mode(dif)[j, d] + 1L, NULL,
           matrix(rowSums(S[, cols, drop = FALSE]), ncol = 1L),
           S[, cols, drop = FALSE])
  }))
}
