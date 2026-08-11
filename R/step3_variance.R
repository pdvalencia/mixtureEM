# ==============================================================================
# Standard errors for step-three covariate models
# ==============================================================================
#
# A three-step latent class analysis relates class membership to covariates in a
# model whose "data" — the posterior class assignments and, under the ML
# adjustment, the classification table built from them — are themselves
# estimates from step one. Two things follow, and the package used to ignore
# both:
#
#   1. The Hessian `m_step.covariate()` stores is the curvature of the
#      *Q function*, the expected complete-data log-likelihood. Under the ML
#      adjustment the marginal step-three log-likelihood is less curved than
#      that (Louis, 1982), so inverting the Q Hessian understates the variance.
#      And with proportional assignment each case contributes K weighted
#      records, so even the correct information matrix has to be sandwiched to
#      account for the duplication (Vermunt, 2010; Bakk, Oberski & Vermunt,
#      2014). Robust standard errors are the standard treatment for exactly
#      this reason.
#
#   2. Holding step one fixed treats estimated quantities as known. The variance
#      of the third-step estimates therefore carries a second component, the
#      variance propagated from step one (Gong & Samaniego, 1981; Oberski &
#      Satorra, 2013; Bakk, Oberski & Vermunt, 2014).
#
# Writing theta1 for the step-one parameters, theta3 for the multinomial-logit
# coefficients, and L3 for the third-step log-likelihood, this file computes
#
#     D3* = D3 + J D1 J' ,      J = (-H3)^{-1} C ,   C = d^2 L3 / d theta3 d theta1'
#
# with H3 the marginal (not Q-function) Hessian of L3, D3 = (-H3)^{-1} M (-H3)^{-1}
# the sandwich whose meat M aggregates the case-level scores (PSU-level under a
# survey design), and D1 the sampling variance of the step-one estimates. This is
# the *first-order* corrected estimator of Bakk, Oberski & Vermunt (2014, eq. 17)
# — the one their simulation recommends, being both better behaved than the
# uncorrected variance and much cheaper than the second-order form of their
# eq. 18, whose extra term vanishes asymptotically (Parke, 1986).
#
# References
#   Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#     assignments to external variables: Standard errors for correct inference.
#     Political Analysis, 22(4), 520-540.
#   Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#     three-step approaches. Political Analysis, 18, 450-469.

#' Standard Errors for Covariate Effects in Three-Step Models
#'
#' @description
#' Relating latent class membership to covariates in two or three steps makes
#' the third step's "data" — the posterior class assignments, and under the ML
#' adjustment the classification table built from them — estimates rather than
#' observations. Ignoring that gives confidence intervals which are too narrow.
#' The \code{se} argument of \code{\link{fit_mixture}} selects among three
#' estimators of the covariate covariance matrix.
#'
#' @details
#' \describe{
#'   \item{\code{"corrected"} (default)}{The first-order corrected estimator of
#'     Bakk et al. (2014, eq. 17),
#'     \eqn{D_3^* = D_3 + J D_1 J'}, where \eqn{D_3} is the step-3 sandwich
#'     below, \eqn{D_1} the sampling variance of the step-1 measurement model,
#'     and \eqn{J = (-H_3)^{-1} \partial^2 L_3 / \partial\theta_3
#'     \partial\theta_1'} the implicit derivative of the step-3 estimates with
#'     respect to the step-1 parameters. This is the estimator their simulation
#'     recommends; the second-order form of their eq. 18 adds a term that
#'     vanishes asymptotically.}
#'   \item{\code{"robust"}}{The step-3 sandwich only,
#'     \eqn{D_3 = (-H_3)^{-1} M (-H_3)^{-1}}, with \eqn{H_3} the marginal
#'     step-3 Hessian and \eqn{M} the outer product of the case-level scores
#'     (PSU-level within strata when a survey design is attached). Proportional
#'     assignment gives each case \eqn{K} weighted records, so this sandwich is
#'     needed even before any step-1 uncertainty is considered.}
#'   \item{\code{"hessian"}}{\eqn{(-H_3)^{-1}} alone. Provided for comparison
#'     with software that reports it; it ignores the record duplication and is
#'     not recommended.}
#' }
#'
#' None of these is the Hessian the package reported previously, which came from
#' the M-step and measured the curvature of the \emph{Q function} rather than of
#' the step-3 log-likelihood, and so was smaller still.
#'
#' How much this matters depends almost entirely on how well separated the
#' classes are. A 250-replication coverage study on the design of Bakk et al.
#' (\code{data-raw/covariate_se_simulation.R} in the package sources) gives, for
#' a nominal 95 percent interval, ranged over six covariate coefficients:
#'
#' \tabular{lllll}{
#'   \strong{design} \tab \strong{entropy R2} \tab \strong{Q Hessian}
#'     \tab \strong{robust} \tab \strong{corrected} \cr
#'   n = 500, rho = .80  \tab .64 \tab .54 - .88 \tab .83 - .98 \tab .92 - .98 \cr
#'   n = 500, rho = .90  \tab .88 \tab .87 - .94 \tab .94 - .98 \tab .94 - .98 \cr
#'   n = 2000, rho = .90 \tab .88 \tab .87 - .95 \tab .94 - .97 \tab .94 - .97
#' }
#'
#' where the first column is the estimator this package reported before the
#' three above existed. With poorly separated classes the correction is the
#' difference between an interval that covers and one that does not, and it bites
#' hardest on the largest effects. Once entropy R-squared reaches about .90 the
#' corrected and robust estimators agree to the third decimal, which matches Bakk
#' et al.'s finding that above n = 2000 with entropy over .90 no correction is
#' needed.
#'
#' A caveat the same study makes visible: at low separation the three-step
#' \emph{point} estimates are themselves biased (Bakk et al., discussion), because
#' step 1 underestimates the classification error. No variance estimator repairs
#' that. A covariate coefficient from a model with entropy below about .60 should
#' be treated with suspicion however wide its interval.
#'
#' \code{summary()}, \code{\link{confint.mixture_model}} and
#' \code{\link{analytical_wald_test}} all print the name of the estimator that
#' produced the numbers they show.
#'
#' @section Scope:
#' The corrected and robust estimators cover a covariate (class-prediction)
#' structural model estimated with \code{n_steps = 2}, or with
#' \code{n_steps = 3} and \code{correction = "none"} or \code{"ML"}. Three cases
#' fall back to the uncorrected Hessian, and say so in the printed output:
#' \code{correction = "BCH"} (whose weights need their own variance treatment,
#' and which is not recommended for covariates in any case);
#' \code{n_steps = 1}, where measurement and structural parameters are estimated
#' jointly and no carry-over correction applies; and a covariate combined with a
#' distal outcome in one nested structural model. The step-1 term additionally
#' requires a measurement model whose parameters this package can put on an
#' unconstrained scale — binary, polytomous, Gaussian, count, mixed, and
#' repeated-measures models qualify; growth models do not, and there the robust
#' sandwich is reported instead.
#'
#' @references
#' Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent class
#' assignments to external variables: Standard errors for correct inference.
#' \emph{Political Analysis}, \emph{22}(4), 520-540. \doi{10.1093/pan/mpu003}
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. \emph{Political Analysis}, \emph{18}(4), 450-469.
#'
#' Gong, G., & Samaniego, F. J. (1981). Pseudo maximum likelihood estimation:
#' Theory and applications. \emph{The Annals of Statistics}, \emph{9}(4), 861-869.
#'
#' @seealso \code{\link{fit_mixture}}, \code{\link{bootstrap_covariates}},
#'   \code{\link{confint.mixture_model}}, \code{\link{analytical_wald_test}}
#'
#' @name covariate_se
NULL

# ------------------------------------------------------------------------------
# 1. The step-three likelihood, its score and its Hessian
# ------------------------------------------------------------------------------
#
# Under the ML adjustment the third-step log-likelihood is
#
#     L3 = sum_i w_i sum_k r_ik log Z_ik ,   Z_ik = sum_j p_j(z_i) C_jk
#
# where r is the matrix of step-one posteriors, C_jk = P(assigned k | true j) is
# the row-normalised classification table, and p_j(z_i) the multinomial logit.
# Writing
#
#     A_ij     = sum_k r_ik C_jk / Z_ik           W_ij = p_ij A_ij
#     B_i(j,l) = sum_k r_ik C_jk C_lk / Z_ik^2
#
# the score and Hessian in the free coefficients (class K anchored at zero) are
#
#     s_i(beta_j) = w_i z_i (W_ij - p_ij)
#     H3[j,l]     = sum_i w_i z_i z_i' [ 1{j=l} (W_ij - p_ij)
#                                        + p_ij p_il (1 - B_i(j,l)) ] .
#
# W sums to one across j, so the score is the ordinary multinomial-logit score
# with W in place of an observed outcome — which is why the EM in `fit_ml()`
# converges to a root of it. With C = I (no adjustment for classification error)
# W collapses to r, B_i(j,l) to 1{j=l} r_ij / p_ij^2, and both expressions
# collapse to the ordinary weighted multinomial logit, which is why the
# unadjusted three-step path shares this code.

.step3_pieces <- function(beta, Zmat, resp1, Cn) {
  P <- softmax_rows(Zmat %*% t(beta))
  if (is.null(Cn)) {
    # Unadjusted: no classification-error correction. Skipping the (numerically
    # awkward) division by P keeps small class probabilities harmless.
    return(list(P = P, W = resp1, Zk = P, Cn = NULL))
  }
  Zk <- pmax(P %*% Cn, 1e-300)
  list(P = P, W = P * ((resp1 / Zk) %*% t(Cn)), Zk = Zk, Cn = Cn)
}

# n x ((K-1)*D) matrix of case-level scores, blocked by class.
.step3_scores <- function(pieces, Zmat, w, K) {
  do.call(cbind, lapply(seq_len(K - 1L), function(j)
    sweep(Zmat, 1, w * (pieces$W[, j] - pieces$P[, j]), "*")))
}

# ((K-1)*D)^2 Hessian of the marginal step-three log-likelihood.
.step3_hessian <- function(pieces, Zmat, w, resp1, K) {
  D <- ncol(Zmat)
  H <- matrix(0, (K - 1L) * D, (K - 1L) * D)
  RZ2 <- if (is.null(pieces$Cn)) NULL else resp1 / pieces$Zk^2

  for (j in seq_len(K - 1L)) {
    for (l in j:(K - 1L)) {
      if (is.null(pieces$Cn)) {
        # B_i(j,l) = 1{j=l} r_ij / p_ij^2, so p_ij p_il (1 - B) is
        # p_ij p_il - 1{j=l} r_ij, and the diagonal term cancels r out entirely.
        coef <- w * (pieces$P[, j] * pieces$P[, l] -
                       if (j == l) pieces$P[, j] else 0)
      } else {
        Bjl  <- as.vector(RZ2 %*% (pieces$Cn[j, ] * pieces$Cn[l, ]))
        coef <- w * ((if (j == l) pieces$W[, j] - pieces$P[, j] else 0) +
                       pieces$P[, j] * pieces$P[, l] * (1 - Bjl))
      }
      blk <- t(Zmat) %*% sweep(Zmat, 1, coef, "*")
      H[((j - 1L) * D + 1L):(j * D), ((l - 1L) * D + 1L):(l * D)] <- blk
      if (l != j)
        H[((l - 1L) * D + 1L):(l * D), ((j - 1L) * D + 1L):(j * D)] <- t(blk)
    }
  }
  H
}

# The classification table P(assigned = k | true = j) under proportional
# assignment, row-normalised, as in `fit_ml()`.
.step3_classification_table <- function(resp1, w) {
  sweep(t(resp1 * w) %*% resp1, 1, colSums(resp1 * w), "/")
}

# ------------------------------------------------------------------------------
# 2. Packing the step-one parameters
# ------------------------------------------------------------------------------
#
# The correction differentiates the third-step score with respect to the
# measurement-model parameters, which therefore need an unconstrained vector
# representation: logits for probabilities, logs for rates and standard
# deviations, means as they are, and class weights as log-ratios against the
# last class. Anything without a representation here returns NULL, and the
# caller falls back to the uncorrected robust variance rather than guessing.
#
# `items` selects a subset of indicators, which is what lets a `blocks` model
# pack a parameter held equal across occasions exactly once.

.step1_n_items <- function(emis) {
  p <- emis$parameters
  m <- p$pis %||% p$means %||% p$rates
  if (is.null(m)) return(NA_integer_)
  if (!is.null(emis$max_val)) ncol(m) %/% emis$max_val else ncol(m)
}

.step1_families <- c("bernoulli", "multinoulli", "poisson",
                     "gaussian_unit", "gaussian_diag")

# The family name with any missing-data suffix removed: the `_nan` variants
# differ only in how the likelihood masks unobserved cells, never in what the
# parameters are, so they pack identically.
.step1_family <- function(emis) sub("_nan$", "", class(emis)[1])

.step1_pack_sub <- function(emis, items = NULL) {
  fam <- .step1_family(emis)
  p   <- emis$parameters
  # Growth and other structured emissions carry parameters this scheme has no
  # unconstrained representation for; they return NULL so the caller can fall
  # back to the uncorrected sandwich instead of guessing.
  if (!fam %in% .step1_families || is.na(.step1_n_items(emis))) return(NULL)
  if (is.null(items)) items <- seq_len(.step1_n_items(emis))
  if (!length(items)) return(numeric(0))
  cols <- .item_param_cols(emis, items)

  switch(fam,
    bernoulli = qlogis(pmin(pmax(as.vector(p$pis[, cols, drop = FALSE]),
                                 1e-10), 1 - 1e-10)),
    poisson   = log(pmax(as.vector(p$rates[, cols, drop = FALSE]), 1e-12)),
    gaussian_unit = as.vector(p$means[, cols, drop = FALSE]),
    gaussian_diag = c(as.vector(p$means[, cols, drop = FALSE]),
                      0.5 * log(pmax(as.vector(p$covariances[, cols, drop = FALSE]),
                                     1e-12))),
    multinoulli = {
      M <- emis$max_val
      unlist(lapply(items, function(j) {
        cj  <- ((j - 1L) * M + 1L):(j * M)
        Pj  <- pmax(p$pis[, cj, drop = FALSE], 1e-12)
        as.vector(log(Pj[, -M, drop = FALSE] / Pj[, M]))
      }), use.names = FALSE)
    },
    NULL)
}

.step1_unpack_sub <- function(emis, par, items = NULL) {
  fam <- .step1_family(emis)
  if (!fam %in% .step1_families || is.na(.step1_n_items(emis))) return(emis)
  if (is.null(items)) items <- seq_len(.step1_n_items(emis))
  if (!length(items)) return(emis)
  cols <- .item_param_cols(emis, items)
  K    <- emis$n_components
  nc   <- length(cols)

  if (fam == "bernoulli") {
    emis$parameters$pis[, cols] <- matrix(plogis(par), K, nc)
  } else if (fam == "poisson") {
    emis$parameters$rates[, cols] <- matrix(exp(par), K, nc)
  } else if (fam == "gaussian_unit") {
    emis$parameters$means[, cols] <- matrix(par, K, nc)
  } else if (fam == "gaussian_diag") {
    emis$parameters$means[, cols] <- matrix(par[seq_len(K * nc)], K, nc)
    emis$parameters$covariances[, cols] <-
      matrix(exp(2 * par[(K * nc + 1L):(2L * K * nc)]), K, nc)
  } else if (fam == "multinoulli") {
    M   <- emis$max_val
    per <- K * (M - 1L)
    for (i in seq_along(items)) {
      j  <- items[i]
      cj <- ((j - 1L) * M + 1L):(j * M)
      L  <- cbind(matrix(par[((i - 1L) * per + 1L):(i * per)], K, M - 1L), 0)
      emis$parameters$pis[, cj] <- softmax_rows(L)
    }
  }
  emis
}

# Which items each block of a `blocks` model contributes freely: block one
# carries every item, later blocks only the ones not held equal across blocks.
.blocks_free_items <- function(mm) {
  lapply(seq_len(mm$n_blocks), function(b)
    if (b == 1L) seq_len(mm$n_items) else
      setdiff(seq_len(mm$n_items), mm$invariant_items))
}

# Pack an emission (measurement model) — NULL when the family is unsupported.
.step1_pack_mm <- function(mm) {
  if (inherits(mm, "blocks")) {
    free  <- .blocks_free_items(mm)
    parts <- lapply(seq_along(mm$models),
                    function(b) .step1_pack_sub(mm$models[[b]], free[[b]]))
    if (any(vapply(parts, is.null, logical(1)))) return(NULL)
    return(unlist(parts, use.names = FALSE))
  }
  if (inherits(mm, "nested")) {
    parts <- lapply(mm$models, .step1_pack_sub)
    if (any(vapply(parts, is.null, logical(1)))) return(NULL)
    return(unlist(parts, use.names = FALSE))
  }
  .step1_pack_sub(mm)
}

.step1_unpack_mm <- function(mm, par) {
  if (inherits(mm, "blocks")) {
    free <- .blocks_free_items(mm)
    at   <- 0L
    for (b in seq_along(mm$models)) {
      n <- length(.step1_pack_sub(mm$models[[b]], free[[b]]))
      mm$models[[b]] <- .step1_unpack_sub(mm$models[[b]],
                                          par[at + seq_len(n)], free[[b]])
      at <- at + n
    }
    # Items held equal across blocks are carried by block one alone.
    if (length(mm$invariant_items) && mm$n_blocks > 1L)
      for (b in 2:mm$n_blocks)
        mm$models[[b]] <- .copy_item_params(mm$models[[b]], mm$models[[1]],
                                            mm$invariant_items)
    return(mm)
  }
  if (inherits(mm, "nested")) {
    at <- 0L
    for (nm in names(mm$models)) {
      n <- length(.step1_pack_sub(mm$models[[nm]]))
      mm$models[[nm]] <- .step1_unpack_sub(mm$models[[nm]], par[at + seq_len(n)])
      at <- at + n
    }
    return(mm)
  }
  .step1_unpack_sub(mm, par)
}

# Full step-one parameter vector: the measurement model followed by the K-1
# free class weights on the log-ratio scale.
.step1_pack <- function(model_state) {
  pm <- .step1_pack_mm(model_state$mm)
  if (is.null(pm)) return(NULL)
  wts <- pmax(model_state$weights, 1e-12)
  c(pm, log(wts[-length(wts)] / wts[length(wts)]))
}

.step1_unpack <- function(model_state, par) {
  K  <- model_state$n_components
  nm <- length(par) - (K - 1L)
  model_state$mm <- .step1_unpack_mm(model_state$mm, par[seq_len(nm)])
  lr <- c(par[nm + seq_len(K - 1L)], 0)
  lr <- lr - max(lr)
  model_state$weights <- exp(lr) / sum(exp(lr))
  model_state
}

# ------------------------------------------------------------------------------
# 3. The corrected variance
# ------------------------------------------------------------------------------

# Sampling variance of the step-one estimates.
#
# `hessian` inverts the numerical observed information of the step-one
# log-likelihood — the textbook estimator. It
# costs O(p^2) likelihood evaluations, so for a large measurement model
# `outer` is used instead: the inverse outer product of the case-level scores
# (BHHH), which needs only the first derivatives that
# have already been computed. Bakk, Oberski & Vermunt (2014, Table 3) cross the
# two choices with everything else and find the difference immaterial.
#
# Returns a list of `V` and the `method` actually used, which is not always the
# one asked for: see the fallback below.
.step1_variance <- function(model_state, X, par, scores, w, method) {
  outer_v <- function() pinv(compute_survey_B(scores, rep(1L, nrow(scores)),
                                              seq_len(nrow(scores))))
  if (method == "outer")
    return(list(V = outer_v(), method = "outer", fallback = FALSE))

  p  <- length(par)
  h  <- .step1_fd_step * pmax(1, abs(par))
  ll <- function(v) sum(w * .step1_ll_case(model_state, X, v))

  H <- matrix(0, p, p)
  for (a in seq_len(p)) {
    for (b in a:p) {
      ea <- numeric(p); ea[a] <- h[a]
      eb <- numeric(p); eb[b] <- h[b]
      H[a, b] <- H[b, a] <-
        (ll(par + ea + eb) - ll(par + ea - eb) -
           ll(par - ea + eb) + ll(par - ea - eb)) / (4 * h[a] * h[b])
    }
  }

  V1 <- pinv(-H)
  V1 <- (V1 + t(V1)) / 2
  ev <- eigen(V1, symmetric = TRUE, only.values = TRUE)$values
  if (!all(is.finite(ev)) || min(ev) < -1e-8 * max(1, max(ev))) {
    # Bakk, Oberski & Vermunt (2014, p. 527) define the corrected variance as a
    # sum of two positive-definite terms, so an indefinite step-one variance
    # breaks that guarantee. Fall back to the outer-product estimator, which is
    # positive semi-definite by construction; their Table 3 finds the two
    # choices immaterial, so the fallback costs nothing statistical.
    #
    # The known trigger is an equality constraint on the measurement model that
    # the packing above does not represent. .step1_pack_sub() exposes every
    # K x J variance cell as a free coordinate, so under `variances_equal` the
    # vector carries K times more variance directions than the model actually
    # has, and the fit is not stationary along the ones the constraint removed.
    # Measured on a four-class continuous model (K = 4, J = 3, n ~ 1000): the
    # step-one gradient is ~0.25 in the means and ~0.16 in the weights but 42.7
    # in the log-sd block, and min eig(-H) = -5.25; refitting the same data with
    # free variances gives min eig(-H) = +7.32. Switching every prior off leaves
    # both unchanged (41.9 and -6.40), so this is the constraint and not, as
    # once supposed, the gap between the penalised mode and the unpenalised
    # likelihood differentiated here.
    return(list(V = outer_v(), method = "outer", fallback = TRUE))
  }
  list(V = V1, method = "hessian", fallback = FALSE)
}

# Case-level step-one log-likelihood at an arbitrary theta1.
.step1_ll_case <- function(model_state, X, par) {
  ms <- .step1_unpack(model_state, par)
  logsumexp(sweep(log_likelihood(ms$mm, X), 2,
                  log(pmax(ms$weights, 1e-300)), "+"), MARGIN = 1)
}

#' @keywords internal
#' @noRd
#
# Variance-covariance matrix of a step-three covariate model.
#
# @param model_state Fitted mixture model whose `mm` is the frozen step-one
#   measurement model and whose `sm` is the covariate model.
# @param X Indicator matrix the measurement model was fitted to.
# @param Zmat Step-three design matrix, intercept column included.
# @param resp1 n x K step-one posterior probabilities.
# @param Cn Row-normalised classification table, or NULL for an unadjusted
#   (`correction = "none"`) third step.
# @param w Case weights.
# @param se One of "corrected", "robust", "hessian".
#
# @return A list with `V` (the K*D square covariance matrix, the anchored class
#   padded with zeros) and `method`, a label naming the estimator used.
.step3_covariate_vcov <- function(model_state, X, Zmat, resp1, Cn, w,
                                  se = "corrected") {
  K <- model_state$n_components
  D <- ncol(Zmat)
  sm <- model_state$sm
  # The free parameters are classes 1..K-1, class K being the anchor .fit_mnl()
  # pins at zero. Every caller runs before sort_model_classes(), which is what
  # keeps the anchor in the last row and lets pad() put the zero block there.
  beta <- sm$parameters$beta

  pad <- function(V) {
    out <- matrix(0, K * D, K * D)
    if ((K - 1L) * D > 0L)
      out[seq_len((K - 1L) * D), seq_len((K - 1L) * D)] <- V
    out
  }
  if (K < 2L || D < 1L)
    return(list(V = matrix(0, K * D, K * D), method = "None (no free parameters)"))

  pieces <- .step3_pieces(beta, Zmat, resp1, Cn)
  H3     <- .step3_hessian(pieces, Zmat, w, resp1, K)
  H3_inv <- pinv(-H3)

  if (identical(se, "hessian"))
    return(list(V = pad(H3_inv), method = "Observed information (step 3 only)"))

  # --- D3: sandwich over case-level (or PSU-level) scores ---------------------
  S      <- .step3_scores(pieces, Zmat, w, K)
  strata <- if (isTRUE(model_state$sm$has_survey_design) &&
                length(model_state$sm$strata) == nrow(Zmat))
    model_state$sm$strata else rep(1L, nrow(Zmat))
  cluster <- if (isTRUE(model_state$sm$has_survey_design) &&
                 length(model_state$sm$cluster) == nrow(Zmat))
    model_state$sm$cluster else seq_len(nrow(Zmat))

  D3 <- H3_inv %*% compute_survey_B(S, strata, cluster) %*% H3_inv
  label_D3 <- if (isTRUE(model_state$sm$has_survey_design))
    "survey-linearized" else "robust"

  if (!identical(se, "corrected"))
    return(list(V = pad(D3),
                method = sprintf("Step-3 sandwich (%s)", label_D3)))

  # --- The step-one uncertainty term ------------------------------------------
  th1 <- .step1_pack(model_state)
  if (is.null(th1) || !length(th1))
    return(list(V = pad(D3),
                method = sprintf(paste("Step-3 sandwich (%s); step-1 correction",
                                       "unavailable for this measurement model"),
                                 label_D3)))

  p1 <- length(th1)
  h1 <- .step1_fd_step * pmax(1, abs(th1))

  # Total step-three score as a function of theta1, at theta3 held fixed. Both
  # the posteriors and the classification table move with theta1.
  s3_of <- function(v) {
    ms <- .step1_unpack(model_state, v)
    ms$sm <- NULL
    r <- exp(e_step(ms, X, NULL)$log_resp)
    Cv <- if (is.null(Cn)) NULL else .step3_classification_table(r, w)
    colSums(.step3_scores(.step3_pieces(beta, Zmat, r, Cv), Zmat, w, K))
  }

  Cmat <- vapply(seq_len(p1), function(m) {
    a <- th1; a[m] <- a[m] + h1[m]
    b <- th1; b[m] <- b[m] - h1[m]
    (s3_of(a) - s3_of(b)) / (2 * h1[m])
  }, numeric((K - 1L) * D))
  dim(Cmat) <- c((K - 1L) * D, p1)

  # Case-level step-one scores, by central differences of the case-level
  # log-likelihood: p1 pairs of evaluations for the whole n x p1 matrix.
  S1 <- vapply(seq_len(p1), function(m) {
    a <- th1; a[m] <- a[m] + h1[m]
    b <- th1; b[m] <- b[m] - h1[m]
    w * (.step1_ll_case(model_state, X, a) -
           .step1_ll_case(model_state, X, b)) / (2 * h1[m])
  }, numeric(nrow(X)))
  dim(S1) <- c(nrow(X), p1)

  d1_method <- if (p1 <= .step1_hessian_max) "hessian" else "outer"
  if (d1_method == "outer")
    message(sprintf(paste("Step-1 variance for the corrected standard errors",
                          "uses the outer-product estimator: the measurement",
                          "model has %d parameters (limit %d for the numerical",
                          "Hessian)."), p1, .step1_hessian_max))
  d1 <- .step1_variance(model_state, X, th1, S1, w, d1_method)
  D1 <- d1$V

  J  <- H3_inv %*% Cmat
  V  <- D3 + J %*% D1 %*% t(J)

  list(V = pad(V),
       method = sprintf("Bakk-Oberski-Vermunt corrected (%s step 3, %s step 1)%s",
                        label_D3, d1$method,
                        if (isTRUE(d1$fallback))
                          "; step-1 Hessian was not positive definite" else ""))
}

# Above this many step-one parameters the numerical Hessian's O(p^2) likelihood
# evaluations stop being worth their cost and the outer-product estimator is
# used instead. 100 parameters is ~20,000 evaluations, a few seconds on a
# typical model; the estimators are asymptotically equivalent either way.
.step1_hessian_max <- 100L

# Relative step for every finite difference taken with respect to the step-one
# parameters. They live on logit/log/identity scales, so the step is scaled by
# max(1, |theta|) to stay meaningful for an item probability pushed far out on
# the logit scale.
#
# 1e-5 is chosen from the middle of the flat part of the accuracy curve rather
# than by taking the smallest step available. On the sleep-quality model the
# corrected standard errors are identical to four decimals from 1e-3 down to
# 1e-6 and then break at 1e-7, where the four-point second difference in
# .step1_variance() divides by h^2 = 1e-14 and loses to cancellation.
.step1_fd_step <- 1e-5

# Design matrix a covariate model builds from Y, kept in one place so that the
# variance code and `m_step.covariate()` cannot drift apart.
.covariate_design <- function(sm, Y) {
  Z <- complete_covariates(as.matrix(Y))
  if (isTRUE(sm$intercept)) cbind(1, Z) else Z
}

# Compute and attach the step-three variance to a covariate sub-model. Shared by
# the ML-adjusted path in `fit_ml()` and the unadjusted path in
# `fit_mixture_internal()`; `Cn = NULL` selects the latter.
.attach_step3_covariate_vcov <- function(model_state, X, Y, resp1, Cn, w,
                                         se = "corrected") {
  if (!inherits(model_state$sm, "covariate")) return(model_state)
  if (identical(se, "none")) return(model_state)

  Zmat <- .covariate_design(model_state$sm, Y)
  res  <- tryCatch(
    .step3_covariate_vcov(model_state, X, Zmat, resp1, Cn, w, se = se),
    error = function(e) {
      warning(sprintf(
        "Step-3 covariate standard errors fell back to the Q-function Hessian: %s",
        conditionMessage(e)), call. = FALSE)
      NULL
    })
  if (is.null(res)) return(model_state)

  model_state$sm$parameters$V_robust <- res$V
  model_state$sm$parameters$V_method <- res$method
  model_state
}
