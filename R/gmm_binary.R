# ==============================================================================
# Growth mixture model for a binary outcome
# ==============================================================================
#
# One binary outcome observed at T occasions. Within class k a case's log-odds
# at occasion t is
#
#     eta_it = Lambda_t alpha_k + Lambda_r,t L_k z_i,      z_i ~ N(0, I),
#
# so the random growth factors are b_i = L_k z_i with covariance
# Psi_k = L_k L_k'. This is the latent class growth model of R/lcga.R plus
# random growth factors, as R/structured_normal.R has them for a continuous
# outcome. There the random effects integrate out in closed form; on the logit
# scale they do not, so the class likelihood
#
#     p(y_i | k) = integral prod_t p(y_it | eta_it) phi(z; 0, I) dz
#
# is evaluated by a product Gauss-Hermite rule in the q standardised random
# effects z (.gauss_hermite(), R/quadrature.R).
#
# The rule is the same for every case: Q nodes per random effect, Q^q in
# all. With every case on the same nodes the linear predictor is one Q x T
# matrix per class and the log-likelihood of every response pattern is two
# matrix products, which is what makes the model affordable. The cost is
# accuracy where the random effects are wide: a variance of 5 on the logit
# scale needs about 30 nodes per dimension for the log-likelihood to sit
# within 0.003 of its limit, and a fit whose variance runs larger still can
# find a maximum of the rule's own error rather than of the likelihood
# (RECORDS.md, "Part 54 E1 (W11)"). fit_gmm() therefore rescores every fit on
# a finer rule and warns when the two disagree.
#
# An adaptive rule is available (`adaptive = TRUE`; Liu & Pierce, 1994;
# Pinheiro & Bates, 1995): for each pattern and class the nodes u are moved
# to z = m + S u, with m and S S' the posterior mean and covariance of z for
# that pattern at the previous iteration, and each weight is multiplied by
# phi(z) |S| / phi(u), which leaves the rule exact for the same polynomials.
# On binary outcomes over five occasions it was no more accurate at 15 nodes
# than the shared rule at 30, and several times slower, since each pattern
# then has nodes of its own; it is kept for the constructor, not offered by
# fit_gmm(). With one node there is nothing to adapt, the rule puts every
# case at z = 0, and the model is exactly the latent class growth model.
#
# The parameters are alpha_k and the lower-triangular L_k, not Psi_k itself.
# With z as the missing data, its density has no parameters, and L_k enters
# the linear predictor exactly as a regression coefficient on Lambda_r z does,
# so each M-step is a logistic regression on the (pattern, node, occasion)
# cells, with the nodes' posterior weights as case weights (Bock & Aitkin,
# 1981). A Psi_k on the edge of the positive semi-definite cone (a variance at
# zero, a correlation of one) is an ordinary point, L_k with a zero on its
# diagonal. An update of Psi from the nodes' posterior second moments cannot
# leave such an edge once it gets there: with the nodes on a line, so is every
# moment computed from them.
#
# The response is binary, so a case is its pattern of answers and of observed
# occasions, and there are at most 3^T distinct ones however large n is. Every
# quantity below is computed per pattern and expanded to cases only at the end.

#' Constructor for binary growth mixture models
#'
#' @description
#' Sets up the emission state for a growth mixture model of a binary outcome:
#' each class follows its own polynomial trajectory on the logit scale, and
#' cases vary about it through normally distributed random growth factors that
#' are integrated out by Gauss-Hermite quadrature.
#'
#' @param n_components Integer. The number of latent classes to estimate.
#' @param design Numeric matrix with one row per occasion and one column per
#'   growth factor, as built by the polynomial in the time scores.
#' @param random_effects Which growth factors vary within a class:
#'   `"intercept_slope"`, `"intercept"` or `"all"`.
#' @param psi `"equal"`, `"variances"` or `"free"`, as in [`fit_gmm()`].
#' @param n_quadrature Number of quadrature nodes per random effect.
#' @param adaptive Logical. Keep the standard normal's nodes for every
#'   response pattern (`FALSE`, the default) or centre them on each pattern's
#'   posterior.
#' @param start_pool Optional list of `n_components` x `ncol(design)`
#'   trajectory matrices that half the random starts are drawn around;
#'   `NULL` draws every start from the pooled trajectory and the observed
#'   response patterns.
#' @param ... Additional arguments, ignored.
#'
#' @return A list object of class \code{c("binary_growth", "emission")}.
#' @export
binary_growth_model <- function(n_components, design,
                                random_effects = "intercept_slope",
                                psi = "equal", n_quadrature = 30,
                                adaptive = FALSE, start_pool = NULL, ...) {
  p <- ncol(design)
  if (!psi %in% c("equal", "variances", "free"))
    stop("Unknown `psi` specification: ", psi, call. = FALSE)
  r_cols <- switch(
    random_effects,
    intercept       = 1L,
    intercept_slope = seq_len(min(2L, p)),
    all             = seq_len(p),
    stop("Unknown `random_effects` specification for a binary outcome: ",
         random_effects, call. = FALSE))
  q <- length(r_cols)
  # Class-specific variances with shared covariances are written below in a
  # parameterisation that exists for two random effects; with three or more
  # there is no simple one, and the two other forms remain.
  if (identical(psi, "variances") && q > 2L)
    stop(paste0("psi = \"variances\" with a binary outcome is available for ",
                "at most two random effects. Use psi = \"equal\" or ",
                "\"free\", or random_effects = \"intercept_slope\"."),
         call. = FALSE)
  gh <- .gauss_hermite(n_quadrature)
  U  <- unname(as.matrix(expand.grid(rep(list(gh$z), q))))
  lw <- unname(rowSums(log(as.matrix(expand.grid(rep(list(gh$w), q))))))
  # The product rule's corners carry products of tail weights: at 30 nodes in
  # two dimensions, 408 of the 900 weigh under 1e-14 and 2e-13 together, so
  # dropping them moves no log-likelihood and halves the work.
  keep <- lw > log(1e-14) | length(lw) == 1L
  U <- U[keep, , drop = FALSE]; lw <- lw[keep]

  state <- list(
    n_components   = n_components,
    design         = design,
    r_cols         = r_cols,
    random_effects = random_effects,
    psi_type       = psi,
    psi_equal      = identical(psi, "equal"),
    n_quadrature   = as.integer(n_quadrature),
    adaptive       = isTRUE(adaptive) && n_quadrature > 1L,
    nodes          = U,
    log_node_w     = lw,
    chol_index     = which(lower.tri(diag(q), diag = TRUE), arr.ind = TRUE),
    # The binomial family supplies the prior of R/lcga.R (.lcga_prior_rows()),
    # which this emission shares: it is a prior on each class's trajectory,
    # placed where a case with no random deviation sits.
    fam            = .wglm_family("binomial"),
    parameters     = list(),
    # Per class, the centre and scale of each pattern's nodes, keyed
    # by the pattern; NULL until the first M-step, and for a pattern the fit
    # never saw, which then gets the standard nodes.
    node_centres   = NULL,
    # The random effects make EM slow on this model, as on the continuous one
    # (R/structured_normal.R); it is accelerated (.em_squarem(), R/em_core.R)
    # and the finish takes it the rest of the way.
    em_max_iter    = 5000L,
    em_tol         = .bg_em_tol,
    em_squarem     = TRUE,
    # Trajectories that seed half the random starts (init_params()).
    start_pool     = start_pool,
    # The occasions the trajectory's edge is placed on (.bg_anchor()).
    edge_map       = .bg_anchor(design)
  )
  class(state) <- c("binary_growth", "emission")
  state
}

# EM is stopped a little before the shared rule (.em_tol_unpolished) and the
# finish takes the winner the rest of the way. Every start runs to that rule:
# the basin of a class whose curve is close to a step climbs slowly and pulls
# ahead late, so a first pass that ranks the starts early -- after 30
# iterations, or at a loose rule of 1e-3 -- set aside the starts that reached
# the maximum on two stress sets (RECORDS.md, "Part 54 E1 (W11)"); stopped at
# 1e-3 itself, the finish can stop short along such a class's ridge.
.bg_em_tol <- list(abs = 1e-5, rel = 1e-12)

# ------------------------------------------------------------------------------
# The trajectory's edge
# ------------------------------------------------------------------------------
#
# A class whose answers follow its curve without exception (every case
# answering 1 early and 0 late, say) has no finite maximum: the likelihood
# keeps rising as its curve steepens towards a step. The curve is therefore
# held inside an edge, and the edge is on the class's log-odds at the
# occasions, not on its coefficients. A class that answers 1 until the last
# occasion needs an intercept of |S| times the time of the switch: a maximum
# with I = 76.8, S = -12.6 over time scores 0-6 switches at t = 6.1, and the
# same slope with I held at 25 switches at t = 2, which cost that solution
# 9.4 units of log-likelihood. On the same solutions the log-odds at the
# first and last occasions never needed to pass the edge by more than
# a few thousandths of a unit (RECORDS.md, "Part 54 E1.2b").
#
# The edge is placed at the p occasions spread evenly over the design's rows
# (first and last for a linear curve, where a straight line takes its
# extremes), so the log-odds there are an invertible linear map of the p
# coefficients, eta = A alpha. Beyond them a curve of higher degree can pass
# the edge between the anchors; the edge's job is to stop the walk towards
# infinity, which bounding the curve at p points does. With repeated time
# scores the anchors may not determine the curve, and the edge goes back on
# the coefficients.
#
# At the edge's value the probability is 1 in double precision
# (plogis(37) == 1), so the objective cannot tell the edge from the limit:
# a curve walked further gains less than the log-likelihood's own rounding.
.bg_edge <- 37

.bg_anchor <- function(design) {
  p <- ncol(design)
  at <- unique(round(seq(1, nrow(design), length.out = p)))
  A <- if (length(at) == p) design[at, , drop = FALSE]
  if (is.null(A) || rcond(A) < 1e-8) A <- diag(1, p)
  list(A = A, Ainv = solve(A))
}

# The trajectories `coefs` (K x p) with each class's log-odds at the anchors
# held inside the edge; a class already inside is returned as it is.
.bg_clamp <- function(ms, coefs) {
  em <- ms$edge_map %||% .bg_anchor(ms$design)
  eta <- coefs %*% t(em$A)
  out <- rowSums(abs(eta) > .bg_edge) > 0
  if (any(out))
    coefs[out, ] <- pmin(pmax(eta[out, , drop = FALSE], -.bg_edge), .bg_edge) %*%
      t(em$Ainv)
  coefs
}

# ------------------------------------------------------------------------------
# The covariance parameters
# ------------------------------------------------------------------------------
#
# The free vector `v` and the K lower-triangular factors it describes, with the
# Jacobian of each factor's lower elements (in .chol_index order) in `v`:
#
#   equal      one factor, shared by every class: v is its lower elements;
#   free       one factor per class, one after the other;
#   variances  q = 2 only: L_k = [a_k, 0; c / a_k, d_k], so that Psi_k has
#              variances a_k^2 and c^2 / a_k^2 + d_k^2 and the covariance c,
#              the same in every class; v = (a_1..a_K, d_1..d_K, c). With one
#              random effect there is no covariance, and this is `free`.
.bg_psi_type <- function(ms) {
  if (identical(ms$psi_type, "variances") && length(ms$r_cols) == 1L) "free"
  else ms$psi_type
}

.bg_n_psi <- function(ms) {
  K <- ms$n_components; nl <- nrow(ms$chol_index)
  switch(.bg_psi_type(ms), equal = nl, free = K * nl, variances = 2L * K + 1L)
}

.bg_psi_map <- function(ms, v) {
  K <- ms$n_components; li <- ms$chol_index; nl <- nrow(li)
  q <- length(ms$r_cols); nv <- length(v)
  fill <- function(l) { L <- matrix(0, q, q); L[li] <- l; L }
  out <- list(L = vector("list", K), J = vector("list", K))
  for (k in seq_len(K)) {
    switch(.bg_psi_type(ms),
      equal = { l <- v; J <- diag(1, nl) },
      free  = { at <- (k - 1L) * nl + seq_len(nl)
                l <- v[at]; J <- matrix(0, nl, nv); J[, at] <- diag(1, nl) },
      variances = {
        a <- v[k]; d <- v[K + k]; cc <- v[nv]
        l <- c(a, cc / a, d)
        J <- matrix(0, 3L, nv)
        J[1L, k] <- 1
        J[2L, k] <- -cc / a^2; J[2L, nv] <- 1 / a
        J[3L, K + k] <- 1
      })
    out$L[[k]] <- fill(l); out$J[[k]] <- J
  }
  out
}

# The free vector that reproduces the given factors (the inverse of the map).
.bg_psi_unmap <- function(ms, L) {
  li <- ms$chol_index
  switch(.bg_psi_type(ms),
    equal = L[[1L]][li],
    free  = unlist(lapply(L, function(M) M[li])),
    variances = c(vapply(L, function(M) M[1L, 1L], 0),
                  vapply(L, function(M) M[2L, 2L], 0),
                  L[[1L]][2L, 1L] * L[[1L]][1L, 1L]))
}

# Write `v` into the state, with the factors and Psi_k = L_k L_k' it implies.
.bg_set_psi <- function(ms, v) {
  mp <- .bg_psi_map(ms, v)
  ms$parameters$psi_par <- v
  ms$parameters$chol    <- mp$L
  ms$parameters$psi     <- lapply(mp$L, tcrossprod)
  ms
}

# ------------------------------------------------------------------------------
# Likelihood
# ------------------------------------------------------------------------------

# Cases grouped by their answers and observed occasions: `idx` maps each case
# to its pattern, `key` names the pattern, `Y0` holds the answers with
# unobserved occasions at 0 and `O` marks the observed ones (both P x T).
#
# Each row is read as a number in base 3 (0, 1, unobserved), which is exact
# for T up to about 33 occasions. The grouping depends on the data alone, and
# every E-step, M-step and finish evaluation asks for it again, so the last
# one is kept and reused when the same data come back.
.bg_pattern_cache <- new.env(parent = emptyenv())
.bg_patterns <- function(X) {
  if (identical(X, .bg_pattern_cache$X)) return(.bg_pattern_cache$pat)
  O   <- !is.na(X)
  key <- as.vector(ifelse(O, X, 2) %*% 3^(seq_len(ncol(X)) - 1L))
  first <- !duplicated(key)
  idx <- match(key, key[first])
  Y0  <- replace(X, !O, 0)[first, , drop = FALSE]
  pat <- list(idx = idx, key = key[first], Y0 = Y0,
              O = O[first, , drop = FALSE] * 1)
  .bg_pattern_cache$X   <- X
  .bg_pattern_cache$pat <- pat
  pat
}

# log(1 + exp(x)) without overflow.
.bg_log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.bg_lse_rows <- function(M) {
  m <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  m + log(rowSums(exp(M - m)))
}

# Class k's nodes for every pattern: `Z`, (P Q) x q with the pattern varying
# fastest, and `lw`, their log weights, P x Q. The standard rule's where the
# state holds no centre for a pattern.
.bg_nodes <- function(ms, k, pat) {
  U <- ms$nodes; P <- nrow(pat$Y0); Q <- nrow(U); q <- ncol(U)
  lw <- matrix(ms$log_node_w, P, Q, byrow = TRUE)
  ad <- ms$node_centres
  at <- if (isTRUE(ms$adaptive) && !is.null(ad)) match(pat$key, ad$key)
        else rep(NA_integer_, P)
  if (all(is.na(at)))
    return(list(Z = U[rep(seq_len(Q), each = P), , drop = FALSE], lw = lw))
  m  <- matrix(0, P, q); Sv <- matrix(as.vector(diag(q)), P, q * q, byrow = TRUE)
  ok <- !is.na(at)
  m[ok, ]  <- ad$m[[k]][at[ok], , drop = FALSE]
  Sv[ok, ] <- ad$S[[k]][at[ok], , drop = FALSE]
  Ur <- U[rep(seq_len(Q), each = P), , drop = FALSE]
  pr <- rep(seq_len(P), Q)
  Z  <- m[pr, , drop = FALSE]
  for (a in seq_len(q)) for (b in seq_len(a))
    Z[, a] <- Z[, a] + Sv[pr, (b - 1L) * q + a] * Ur[, b]
  logdet <- rowSums(log(abs(Sv[, (seq_len(q) - 1L) * q + seq_len(q), drop = FALSE])))
  lw <- lw + matrix(-0.5 * rowSums(Z^2) + 0.5 * rowSums(Ur^2), P, Q) + logdet
  list(Z = Z, lw = lw)
}

# Linear predictor of class k at its nodes, (P Q) x T.
.bg_eta_at <- function(ms, k, Z) {
  Lr <- ms$design[, ms$r_cols, drop = FALSE]
  B  <- Z %*% t(ms$parameters$chol[[k]])
  sweep(B %*% t(Lr), 2L, as.vector(ms$design %*% ms$parameters$coefs[k, ]), "+")
}

# Node-wise log of weight times likelihood, P x Q.
.bg_node_lp <- function(ms, E, pat, lw) {
  P  <- nrow(pat$Y0); pr <- rep(seq_len(P), ncol(lw))
  matrix(rowSums(pat$Y0[pr, , drop = FALSE] * E -
                   pat$O[pr, , drop = FALSE] * .bg_log1pexp(E)), P) + lw
}

# Whether class k's nodes are the standard ones for every pattern, in which
# case every pattern shares one Q x T linear predictor and the log-likelihood
# of all of them is two matrix products (.bg_class_lp()).
.bg_fixed_nodes <- function(ms, pat)
  !isTRUE(ms$adaptive) || is.null(ms$node_centres) ||
    !any(pat$key %in% ms$node_centres$key)

# Class k's nodes and node-wise log of weight times likelihood (P x Q), by
# the shared-node route where the nodes are the standard ones: `Z` is then
# Q x q and `E` Q x T, and otherwise (P Q) x q and (P Q) x T.
.bg_class_lp <- function(ms, k, pat) {
  if (.bg_fixed_nodes(ms, pat)) {
    E  <- .bg_eta_at(ms, k, ms$nodes)
    lp <- pat$Y0 %*% t(E) - pat$O %*% t(.bg_log1pexp(E))
    return(list(Z = ms$nodes, E = E, fixed = TRUE,
                lp = sweep(lp, 2L, ms$log_node_w, "+")))
  }
  nd <- .bg_nodes(ms, k, pat)
  E  <- .bg_eta_at(ms, k, nd$Z)
  list(Z = nd$Z, E = E, fixed = FALSE, lp = .bg_node_lp(ms, E, pat, nd$lw))
}

# Pattern-by-class log-likelihood (P x K).
.bg_pattern_ll <- function(ms, pat)
  vapply(seq_len(ms$n_components), function(k)
    .bg_lse_rows(.bg_class_lp(ms, k, pat)$lp), numeric(nrow(pat$Y0)))

#' @exportS3Method
log_likelihood.binary_growth <- function(model_state, X, ...) {
  pat <- .bg_patterns(X)
  .bg_pattern_ll(model_state, pat)[pat$idx, , drop = FALSE]
}

#' @exportS3Method
n_parameters.binary_growth <- function(model_state, ...)
  length(model_state$parameters$coefs) + .bg_n_psi(model_state)

# Starting values. A class whose curve is close to a step (cases that answer 1
# until some occasion and 0 after it, or always 1 in a sample that mostly
# answers 0) sits far from the pooled trajectory on the logit scale, and its
# basin is narrow: perturbing the pooled curve, as init_params.lcga() does,
# almost never lands in it, and with two such classes in one solution, never
# (RECORDS.md, "Part 54 E1.2"). So half the starts take the latent class
# growth model's trajectories at the same number of classes (`start_pool`,
# .bg_start_pool()), which hold such classes when the data have them, jittered
# by 0.5 on the logit scale; the other half perturb the pooled trajectory and
# then, with probability one half per class, replace the class by a curve
# through one observed response pattern, chosen uniformly among the distinct
# patterns so that rare ones are drawn as often as common ones: +-A on the
# logit scale where the pattern answers 1 or 0, A uniform on (2, 10), fitted by
# least squares on its observed occasions. Psi starts diagonal, each random
# effect moving the log-odds by about one unit at the occasion where it weighs
# most. A start has no centres yet, and uses the standard nodes.
#' @exportS3Method
init_params.binary_growth <- function(model_state, X, resp, random_state = NULL,
                                      ...) {
  if (!is.null(random_state)) set.seed(random_state)
  K   <- model_state$n_components
  D   <- model_state$design
  pool <- model_state$start_pool
  if (length(pool) && stats::runif(1) < 0.5) {
    C0 <- pool[[sample.int(length(pool), 1L)]]
    coefs <- C0 + matrix(stats::rnorm(length(C0), 0, 0.5), nrow(C0))
  } else {
    lng <- .lcga_long(model_state, X)
    pooled <- .wglm_fit(lng$D, lng$y, as.numeric(lng$observed), model_state$fam)
    eta0 <- as.vector(D %*% pooled$coefficients)
    coefs <- matrix(0, K, ncol(D))
    for (k in seq_len(K))
      coefs[k, ] <- qr.solve(D, eta0 + stats::rnorm(length(eta0)))
    pat <- .bg_patterns(X)
    for (k in seq_len(K)) {
      if (stats::runif(1) >= 0.5) next
      i <- sample.int(nrow(pat$Y0), 1L)
      o <- pat$O[i, ] > 0
      if (sum(o) < ncol(D)) next
      A <- stats::runif(1, 2, 10)
      coefs[k, ] <- qr.solve(D[o, , drop = FALSE], ifelse(pat$Y0[i, o] > 0, A, -A))
    }
  }
  coefs <- .bg_clamp(model_state, coefs)
  Lr <- D[, model_state$r_cols, drop = FALSE]
  L0 <- diag(1 / sqrt(pmax(apply(Lr^2, 2L, max), 1e-8)), ncol(Lr))
  model_state$parameters$coefs <- coefs
  model_state$node_centres <- NULL
  .bg_set_psi(model_state, .bg_psi_unmap(model_state, rep(list(L0), K)))
}

# ------------------------------------------------------------------------------
# M-step
# ------------------------------------------------------------------------------

# The cells' summaries one class contributes, for weights `Wy` (posterior
# weight times answer) and `Wo` (times observed), (P Q) x T at the nodes `Z`:
# its value, and the gradient and curvature of the expected complete-data
# log-likelihood in (alpha_k, L_k's lower elements). A column of the cell
# design is Lambda_t (alpha) or Lambda_r,t[a] z[b] (element (a, b) of L), so
# every sum it needs is a time-by-time sum of the residual or the curvature
# times 1, z_b or z_b z_d, taken once over the cells.
.bg_cell_derivs <- function(ms, k, Z, Wy, Wo, deriv = TRUE) {
  E   <- .bg_eta_at(ms, k, Z)
  val <- sum(Wy * E - Wo * .bg_log1pexp(E))
  if (!deriv) return(list(value = val))
  D  <- ms$design; Lr <- D[, ms$r_cols, drop = FALSE]
  li <- ms$chol_index; p <- ncol(D); nl <- nrow(li)
  mu <- stats::plogis(E)
  R  <- Wy - Wo * mu
  C  <- Wo * mu * (1 - mu)
  # Column j: its time profile and which z it multiplies (0 for none).
  prof <- cbind(D, Lr[, li[, 1L], drop = FALSE])
  zix  <- c(rep(0L, p), li[, 2L])
  # The time-by-time sums, all at once: Z1 = [1, z] against the residual, and
  # every product of two of its columns against the curvature.
  q  <- ncol(Z); nj <- p + nl
  Z1 <- cbind(1, Z)
  pq <- which(upper.tri(diag(q + 1L), diag = TRUE), arr.ind = TRUE)
  M1 <- crossprod(Z1, R)                                     # (q + 1) x T
  M2 <- crossprod(Z1[, pq[, 1L], drop = FALSE] * Z1[, pq[, 2L], drop = FALSE], C)
  slot <- matrix(0L, q + 1L, q + 1L)
  slot[pq] <- seq_len(nrow(pq)); slot[pq[, 2:1, drop = FALSE]] <- seq_len(nrow(pq))
  zi <- zix + 1L
  g <- rowSums(t(prof) * M1[zi, , drop = FALSE])
  H <- matrix(0, nj, nj)
  for (j in seq_len(nj)) for (l in seq_len(j))
    H[j, l] <- H[l, j] <- sum(prof[, j] * prof[, l] * M2[slot[zi[j], zi[l]], ])
  list(value = val, g = g, H = H)
}

# The objective over all classes in (alpha, v), with each class's cells
# (`cells[[k]]`: nodes Z and weights Wy, Wo) and the prior's pseudo-rows `pr`
# (.lcga_prior_rows()) on each class's trajectory; with `deriv`, the gradient
# and curvature, exact in (alpha, L) and carried to v through the map's
# Jacobian (Gauss-Newton where the map is not linear).
.bg_q <- function(ms, cells, pr, deriv = FALSE) {
  K  <- ms$n_components; p <- ncol(ms$design)
  nl <- nrow(ms$chol_index); nv <- length(ms$parameters$psi_par)
  mp <- if (deriv) .bg_psi_map(ms, ms$parameters$psi_par)
  val <- 0
  if (deriv) { g <- numeric(K * p + nv); H <- matrix(0, K * p + nv, K * p + nv) }
  for (k in seq_len(K)) {
    a  <- ms$parameters$coefs[k, ]
    ck <- cells[[k]]
    d  <- .bg_cell_derivs(ms, k, ck$Z, ck$Wy, ck$Wo, deriv)
    val <- val + d$value
    ep <- if (length(pr$y)) as.vector(pr$D %*% a)
    if (length(pr$y)) val <- val + sum(pr$w * (pr$y * ep - .bg_log1pexp(ep)))
    if (!deriv) next
    gk <- d$g; Hk <- d$H
    ia <- seq_len(p); il <- p + seq_len(nl)
    if (length(pr$y)) {
      m <- stats::plogis(ep)
      gk[ia] <- gk[ia] + crossprod(pr$D, pr$w * (pr$y - m))
      Hk[ia, ia] <- Hk[ia, ia] + crossprod(pr$D, pr$w * m * (1 - m) * pr$D)
    }
    J  <- mp$J[[k]]
    # coefs is packed by column (as.vector()), so class k's coefficient j is
    # coordinate (j - 1) K + k.
    ak <- (ia - 1L) * K + k; iv <- K * p + seq_len(nv)
    g[ak] <- gk[ia]
    g[iv] <- g[iv] + crossprod(J, gk[il])
    H[ak, ak] <- Hk[ia, ia]
    H[ak, iv] <- Hk[ia, il, drop = FALSE] %*% J
    H[iv, ak] <- t(H[ak, iv])
    H[iv, iv] <- H[iv, iv] + crossprod(J, Hk[il, il, drop = FALSE] %*% J)
  }
  if (!deriv) return(val)
  list(value = val, g = g, H = H)
}

# Each pattern's posterior mean and covariance of z in one class, from its
# node weights `Hp` (P x Q) at nodes `Z`: the centres the next E-step's rule
# is built on, with S the lower Cholesky factor of the covariance, flattened
# by column. The covariance is that of z given the answers, whose prior is the
# identity, so it is positive definite away from rounding; a pattern where it
# is not keeps the standard scale.
.bg_centres <- function(Hp, Z, q) {
  P <- nrow(Hp); Q <- ncol(Hp); pr <- rep(seq_len(P), Q)
  h <- as.vector(Hp)
  m <- vapply(seq_len(q), function(a) rowsum(h * Z[, a], pr, reorder = TRUE)[, 1L],
              numeric(P))
  m <- matrix(m, P, q)
  S <- matrix(as.vector(diag(q)), P, q * q, byrow = TRUE)
  for (i in seq_len(P)) {
    w  <- h[pr == i]; zi <- Z[pr == i, , drop = FALSE]
    Ci <- crossprod(zi, w * zi) - tcrossprod(m[i, ])
    ch <- tryCatch(t(chol((Ci + t(Ci)) / 2)), error = function(e) NULL)
    if (!is.null(ch) && all(is.finite(ch)) && min(diag(ch)) > 1e-6)
      S[i, ] <- as.vector(ch)
  }
  list(m = m, S = S)
}

#' @exportS3Method
m_step.binary_growth <- function(model_state, X, resp, weights = NULL, ...) {
  if (!is.null(weights)) resp <- sweep(resp, 1, weights, "*")
  ms  <- model_state
  K   <- ms$n_components; p <- ncol(ms$design); q <- length(ms$r_cols)
  # The edge holds whatever brought the state here, an accelerated step
  # (.em_squarem()) included, not only this M-step's own moves.
  ms$parameters$coefs <- .bg_clamp(ms, ms$parameters$coefs)
  pat <- .bg_patterns(X)
  P   <- nrow(pat$Y0)
  R   <- rowsum(resp, pat$idx, reorder = TRUE)   # P x K, patterns in order
  pr  <- .lcga_prior_rows(ms, .lcga_marginal(X, weights))

  # The posterior over each class's nodes, at the parameters (and the centres)
  # the E-step used; the nodes stay where they are for the whole M-step.
  cells <- vector("list", K)
  adapt <- list(key = pat$key, m = vector("list", K), S = vector("list", K))
  for (k in seq_len(K)) {
    cl <- .bg_class_lp(ms, k, pat)
    Hp <- exp(cl$lp - .bg_lse_rows(cl$lp))                # P x Q
    # Shared nodes: each node's cells summed over the patterns first.
    cells[[k]] <- if (cl$fixed) {
      G <- Hp * R[, k]
      list(Z = cl$Z, Wy = crossprod(G, pat$Y0), Wo = crossprod(G, pat$O))
    } else {
      G   <- as.vector(Hp * R[, k])
      pr_ <- rep(seq_len(P), ncol(Hp))
      list(Z = cl$Z, Wy = G * pat$Y0[pr_, , drop = FALSE],
           Wo = G * pat$O[pr_, , drop = FALSE])
    }
    if (isTRUE(ms$adaptive)) {
      Zc <- if (cl$fixed) ms$nodes[rep(seq_len(nrow(ms$nodes)), each = P), , drop = FALSE]
            else cl$Z
      ce <- .bg_centres(Hp, Zc, q)
      adapt$m[[k]] <- ce$m; adapt$S[[k]] <- ce$S
    }
  }

  # One Newton step with step halving, a generalised EM step: the objective
  # is concave in (alpha, L), and under `variances`, where the map to L bends,
  # halving keeps the step an ascent. Iterating the M-step to its maximum gave
  # the same fits at a third more cost: the acceleration and the finish do the
  # rest (RECORDS.md, "Part 54 E1 (W11)").
  get <- function(s) c(as.vector(s$parameters$coefs), s$parameters$psi_par)
  # The trajectory is held inside its edge (.bg_clamp()): a class whose
  # curve is close to a step would otherwise be walked off towards infinity
  # by EM itself.
  put <- function(s, th) {
    s$parameters$coefs[] <- .bg_clamp(s, matrix(th[seq_len(K * p)], K))
    .bg_set_psi(s, th[-seq_len(K * p)])
  }
  cur <- .bg_q(ms, cells, pr, deriv = TRUE)
  for (it in seq_len(1L)) {
    step <- tryCatch(solve(cur$H + diag(1e-8, nrow(cur$H)), cur$g),
                     error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
    th <- get(ms); s <- 1; moved <- FALSE
    repeat {
      cand <- put(ms, th + s * step)
      f1 <- .bg_q(cand, cells, pr)
      if (is.finite(f1) && f1 >= cur$value - 1e-12) { moved <- TRUE; break }
      s <- s / 2
      if (s < 1e-6) break
    }
    if (!moved) break
    ms  <- cand
    if (max(abs(s * step)) < 1e-8) break
    cur <- .bg_q(ms, cells, pr, deriv = TRUE)
  }
  # The next E-step builds its rule on the posterior just used.
  if (isTRUE(ms$adaptive)) ms$node_centres <- adapt
  ms
}

# The class trajectory on the probability scale, K x T: at the class mean of
# the growth factors (the default) and averaged over the random effects
# (`marginal`, the proportion of the class expected to answer 1), on the
# standard nodes, since this average belongs to no response pattern.
.bg_fitted <- function(mm, marginal = FALSE) {
  K <- mm$n_components
  if (!marginal)
    return(t(stats::plogis(mm$design %*% t(mm$parameters$coefs))))
  w <- exp(mm$log_node_w)
  t(vapply(seq_len(K), function(k)
    as.vector(crossprod(stats::plogis(.bg_eta_at(mm, k, mm$nodes)), w)),
    numeric(nrow(mm$design))))
}

# ------------------------------------------------------------------------------
# Packing for the finish and the standard errors (R/qn_pack.R)
# ------------------------------------------------------------------------------

# Each class's log-odds at the edge's anchors (.bg_anchor()), then the
# covariance's free vector `v`. The reported growth means, and Psi, are
# carried to by the delta method. The nodes' centres are not parameters: they
# stay where EM left them, which makes the finish's objective one fixed rule.
#
# The trajectory is packed by its log-odds rather than its coefficients so
# that the finish holds a class at its wall where EM holds it at its edge:
# a class whose curve is close to a step is held at the occasion where it is
# certain, as a response probability run to 0 or 1 is, and moves freely at
# the occasion where it is not, rather than being walked towards infinity
# with a random-effect variance in tow.
.qn_pack_binary_growth <- function(ms) {
  em <- ms$edge_map %||% .bg_anchor(ms$design)
  .qn_join(list(.qn_part(ms$parameters$coefs %*% t(em$A), "logit"),
                .qn_part(ms$parameters$psi_par, "free")))
}

# Each case's score in that packed order, for the finish's analytic route
# (R/qn_scores.R): by Fisher's identity, sum_k r_ik times the node-posterior
# average of "observed minus expected" times the cell design, per pattern, with
# L_k's columns carried to v through the map's Jacobian.
.qn_bg_scores <- function(ms, X, R) {
  K   <- ms$n_components; p <- ncol(ms$design)
  pat <- .bg_patterns(X); P <- nrow(pat$Y0)
  D   <- ms$design; Lr <- D[, ms$r_cols, drop = FALSE]
  li  <- ms$chol_index; nl <- nrow(li); nv <- length(ms$parameters$psi_par)
  prof <- cbind(D, Lr[, li[, 1L], drop = FALSE])
  zix  <- c(rep(0L, p), li[, 2L])
  mp  <- .bg_psi_map(ms, ms$parameters$psi_par)
  em  <- ms$edge_map %||% .bg_anchor(D)
  S   <- matrix(0, nrow(X), K * p + nv)
  for (k in seq_len(K)) {
    cl <- .bg_class_lp(ms, k, pat)
    Hm <- exp(cl$lp - .bg_lse_rows(cl$lp))                # P x Q
    mu <- stats::plogis(cl$E)
    Sp <- if (cl$fixed) {
      # sum_q H_pq z_q (y_pt - mu_qt), observed cells only.
      vapply(seq_len(p + nl), function(j) {
        z <- if (zix[j] == 0L) rep(1, nrow(cl$Z)) else cl$Z[, zix[j]]
        as.vector((pat$O * (pat$Y0 * as.vector(Hm %*% z) -
                              Hm %*% (z * mu))) %*% prof[, j])
      }, numeric(P))
    } else {
      pr <- rep(seq_len(P), ncol(Hm))
      Rc <- as.vector(Hm) * pat$O[pr, , drop = FALSE] *
        (pat$Y0[pr, , drop = FALSE] - mu)                 # (P Q) x T
      vapply(seq_len(p + nl), function(j) {
        z <- if (zix[j] == 0L) 1 else cl$Z[, zix[j]]
        as.vector(rowsum((z * Rc) %*% prof[, j], pr, reorder = TRUE))
      }, numeric(P))
    }
    # The trajectory's columns carried to the anchors' log-odds, alpha = A^-1 eta.
    Sp <- matrix(Sp, P)
    Sp[, seq_len(p)] <- Sp[, seq_len(p), drop = FALSE] %*% em$Ainv
    Sk <- Sp[pat$idx, , drop = FALSE] * R[, k]
    S[, (seq_len(p) - 1L) * K + k] <- Sk[, seq_len(p)]
    S[, K * p + seq_len(nv)] <- S[, K * p + seq_len(nv)] +
      Sk[, p + seq_len(nl), drop = FALSE] %*% mp$J[[k]]
  }
  S
}

.qn_unpack_binary_growth <- function(ms, par) {
  n  <- length(ms$parameters$coefs)
  em <- ms$edge_map %||% .bg_anchor(ms$design)
  ms$parameters$coefs[] <- matrix(par[seq_len(n)], ms$n_components) %*% t(em$Ainv)
  .bg_set_psi(ms, par[-seq_len(n)])
}

# ------------------------------------------------------------------------------
# fit_gmm(family = "binomial")
# ------------------------------------------------------------------------------

.fit_gmm_binary <- function(prep, n_classes, degree, random_effects, psi,
                            time_scores, n_quadrature, predictors, ...) {
  n_quadrature <- suppressWarnings(as.integer(n_quadrature))
  if (length(n_quadrature) != 1L || is.na(n_quadrature) || n_quadrature < 1L)
    stop("`n_quadrature` must be a single whole number of at least 1.",
         call. = FALSE)
  .lcga_check_outcome(prep$X, "binomial")
  design <- .lcga_design(time_scores, degree)

  pool <- .bg_start_pool(prep$X, n_classes, design, ...)

  fit <- fit_mixture(
    indicators     = prep$X,
    n_classes      = n_classes,
    measurement    = "binary_growth",
    predictors     = predictors,
    design         = design,
    random_effects = random_effects,
    psi            = psi,
    n_quadrature   = n_quadrature,
    start_pool     = pool,
    ...
  )

  mm <- fit$mm
  fit$growth <- list(
    model          = "gmm",
    family         = "binomial",
    degree         = degree,
    random_effects = random_effects,
    time_scores    = time_scores,
    time_labels    = prep$time_labels,
    design         = design,
    factor_names   = colnames(design),
    means          = mm$parameters$coefs,
    covariate_equal = FALSE,
    fitted         = .bg_fitted(mm, marginal = TRUE),
    fitted_at_mean = .bg_fitted(mm),
    psi            = mm$parameters$psi[seq_len(if (mm$psi_equal) 1L else
                                                 n_classes)],
    psi_equal      = mm$psi_equal,
    psi_type       = mm$psi_type,
    random_names   = colnames(design)[mm$r_cols],
    n_quadrature   = n_quadrature,
    wave_missing   = prep$wave_missing
  )
  fit$growth$quadrature_check <- .bg_quadrature_check(fit)
  qc <- fit$growth$quadrature_check
  if (!is.null(qc) && is.finite(qc$difference) &&
      abs(qc$difference) > .bg_quadrature_tol)
    warning(sprintf(paste0(
      "The log-likelihood moves by %.3f when the random effects are ",
      "integrated with %d nodes instead of %d, so at this solution the ",
      "integral is not accurate enough to compare models on. A random ",
      "effect with a variance of %.1f on the logit scale is wide for the ",
      "rule. Refit with n_quadrature = %d or more."),
      qc$difference, qc$nodes, mm$n_quadrature, qc$max_variance,
      qc$nodes), call. = FALSE)
  .check_replication(fit)
  class(fit) <- c("gmm", class(fit))
  fit
}

# The latent class growth model at the same number of classes, whose
# trajectories seed half the starts (init_params.binary_growth()): fitted to
# the same cases with the same weights, priors and seed, on at most 20 starts.
# A fit that continues another (`refine_from`, `start_from`) draws no starts,
# and a search of one start gets none.
.bg_start_pool <- function(X, n_classes, design, ...) {
  a <- list(...)
  n_init <- a$n_init %||% 20L
  if (!is.null(a$refine_from) || !is.null(a$start_from) || n_init < 2L)
    return(NULL)
  a <- a[intersect(names(a), c("random_state", "n_cores", "bayes_constants",
                               "smoothing", "weights", "weight_type"))]
  a$n_init <- min(as.integer(n_init), 20L)
  lc <- tryCatch(suppressWarnings(suppressMessages(do.call(fit_mixture, c(list(
    indicators = X, n_classes = n_classes, measurement = "lcga",
    design = design, family = "binomial"), a)))), error = function(e) NULL)
  if (is.null(lc)) NULL else list(lc$mm$parameters$coefs)
}

# The fitted log-likelihood rescored with 10 more nodes per random effect, at
# the same parameters. The difference is what the rule's error is worth at
# this solution; NULL for one node, which is a different model rather than a
# rough integral.
.bg_quadrature_tol <- 0.01
.bg_quadrature_check <- function(fit) {
  mm <- fit$mm
  if (mm$n_quadrature <= 1L) return(NULL)
  Q2 <- mm$n_quadrature + 10L
  fine <- binary_growth_model(mm$n_components, mm$design, mm$random_effects,
                              mm$psi_type, n_quadrature = Q2)
  fine$parameters <- mm$parameters
  ll <- function(m) {
    L <- sweep(log_likelihood(m, fit$data), 2L,
               log(pmax(fit$weights, 1e-300)), "+")
    sum(fit$sample_weights * logsumexp(L, MARGIN = 1))
  }
  # Only the measurement part is rescored, so a one-step model with an
  # outcome or predictors compares like with like.
  list(nodes = Q2, difference = ll(fine) - ll(mm),
       max_variance = max(vapply(mm$parameters$psi, function(P) max(diag(P)), 0)))
}
