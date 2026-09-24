# ==============================================================================
# Direct covariate effects on binary indicators of a cross-sectional model
# ==============================================================================
#
# The cross-sectional counterpart of R/lta_dif.R. There, `predictors_items`
# lets a covariate act on an item inside each latent status of an LTA; here it
# does the same inside each class of an LCA fitted by fit_mixture():
#
#   logit P(Y_ij = 1 | X = k, z_i) = theta_kj + z_i' delta_j        (uniform)
#                                  = theta_kj + z_i' delta_kj       (by class)
#
# The user names every covariate-item pair and whether its slope differs by
# class; nothing is searched (DECISIONS, "Part 50 covariate selection"). An
# item named in no pair keeps the plain Bernoulli emission's closed-form update,
# so a model with no pairs at all is the unconditional fit exactly.
#
# The two consequences R/lta_dif.R's header states carry over unchanged.
#
# The covariate enters through a small number of distinct PATTERNS. The E-step
# builds one K x J probability table per pattern and fills the rows of that
# pattern from it, so the row-work is the plain emission's; the cost is the
# loop. A many-valued covariate is refused by the same cap,
# `options(mixtureEM.dif_max_patterns = )`.
#
# `parameters$pis` holds the probabilities of a case whose item covariates are
# all zero, and the categorical Bayes prior is applied there, at z = 0. Every
# reader of `pis` -- the log-prior in R/em_core.R, class sorting, profile plots,
# the measurement summaries -- therefore works unchanged and reports the
# z = 0 case, which is what the print method says.
#
# `parameters$dif_slopes` is K x (J * D), one K x J block per covariate, with
# the rows of a uniform slope identical and a pair that was not named held at
# zero. It is a K-row matrix so that the generic class permutations move it
# with `pis`.

# ------------------------------------------------------------------------------
# Construction
# ------------------------------------------------------------------------------

# Turn a plain binary emission into a DIF one. `dif` is what fit_mixture()
# builds: `Z` (n x D, the union of the named covariates), `free` (J x D logical,
# which item-covariate slopes exist) and `by_class` (length J logical).
.dif_emission <- function(mm, dif) {
  if (!class(mm)[1] %in% c("bernoulli", "bernoulli_nan"))
    stop("`predictors_items` is available for binary indicators only.",
         call. = FALSE)
  Z <- as.matrix(dif$Z)
  if (anyNA(Z))
    stop("`predictors_items` covariates may not be missing: the item ",
         "probabilities of a case depend on them, and there is nothing to ",
         "integrate a missing value over.", call. = FALSE)

  codes <- lapply(seq_len(ncol(Z)), function(d) match(Z[, d], unique(Z[, d])))
  key   <- do.call(paste, c(codes, sep = "\r"))
  u     <- !duplicated(key)
  pat   <- match(key, key[u])
  P     <- sum(u)
  cap   <- getOption("mixtureEM.dif_max_patterns", 16L)
  if (P > cap)
    stop(sprintf(paste0(
      "`predictors_items` produced %d distinct covariate patterns, above the ",
      "limit of %d. Each pattern costs one extra probability table in every ",
      "E-step and M-step. Bin or dichotomise the covariate, or raise ",
      "`options(mixtureEM.dif_max_patterns = )` deliberately."), P, cap),
      call. = FALSE)

  mm$dif <- list(Z = Z, Zu = Z[u, , drop = FALSE], pat = pat, P = P,
                 rows = lapply(seq_len(P), function(p) which(pat == p)),
                 free = dif$free, by_class = dif$by_class,
                 covariates = colnames(Z))
  class(mm) <- c("bernoulli_dif", "emission")
  mm
}

# Resolve fit_mixture()'s `predictors_items` (a named list, item -> one-sided
# formula or column names) and `predictors_items_by_class` (item names) into
# the `dif` list .dif_emission() takes. Factors are dummy-coded the way
# `predictors` is, by prepare_covariates().
.dif_spec <- function(spec, by_class, data, X, measurement) {
  if (!is.character(measurement) || length(measurement) != 1L ||
      !measurement %in% .binary_families)
    stop("`predictors_items` is available for binary indicators only ",
         "(`measurement = \"binary\"`).", call. = FALSE)
  items <- colnames(X)
  if (!is.list(spec) || is.null(names(spec)) || any(names(spec) == ""))
    stop("`predictors_items` must be a named list: item name = one-sided ",
         "formula or column names, e.g. list(item3 = ~ female).", call. = FALSE)
  bad <- setdiff(c(names(spec), by_class %||% character(0)), items)
  if (length(bad))
    stop(sprintf("`predictors_items` names %s not among the indicators: %s.",
                 if (length(bad) > 1L) "items" else "an item",
                 paste(bad, collapse = ", ")), call. = FALSE)
  bad <- setdiff(by_class %||% character(0), names(spec))
  if (length(bad))
    stop(sprintf(paste0("`predictors_items_by_class` names %s with no entry ",
                        "in `predictors_items`: %s."),
                 if (length(bad) > 1L) "items" else "an item",
                 paste(bad, collapse = ", ")), call. = FALSE)
  if (is.null(data))
    stop("`predictors_items` is resolved against `data`; supply it.",
         call. = FALSE)

  per_item <- lapply(names(spec), function(it)
    as.matrix(prepare_covariates(
      .columns_from_data(spec[[it]], data, "predictors_items"))))
  if (any(vapply(per_item, nrow, 1L) != nrow(X)))
    stop("`data` must have one row per case of the indicators.", call. = FALSE)
  covs <- unique(unlist(lapply(per_item, colnames)))
  Z    <- matrix(0, nrow(X), length(covs), dimnames = list(NULL, covs))
  free <- matrix(FALSE, length(items), length(covs),
                 dimnames = list(items, covs))
  for (i in seq_along(per_item)) {
    Z[, colnames(per_item[[i]])] <- per_item[[i]]
    free[names(spec)[i], colnames(per_item[[i]])] <- TRUE
  }
  list(Z = Z, free = free, by_class = items %in% (by_class %||% character(0)))
}

# Column of `dif_slopes` holding item j's slope on covariate d.
.dif_col <- function(J, j, d) (d - 1L) * J + j

# The K x J shift every item's logit takes for covariate pattern p.
.dif_shift <- function(mm, p) {
  S  <- mm$parameters$dif_slopes
  J  <- ncol(mm$parameters$pis)
  zu <- mm$dif$Zu
  out <- matrix(0, nrow(S), J)
  for (d in seq_len(ncol(zu)))
    if (zu[p, d] != 0)
      out <- out + S[, (d - 1L) * J + seq_len(J), drop = FALSE] * zu[p, d]
  out
}

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------

#' @exportS3Method
init_params.bernoulli_dif <- function(model_state, X, resp, random_state = NULL, ...) {
  model_state <- init_params.bernoulli(model_state, X, resp, random_state)
  model_state$parameters$dif_slopes <-
    matrix(0, model_state$n_components, ncol(X) * ncol(model_state$dif$Z))
  model_state
}

#' @exportS3Method
log_likelihood.bernoulli_dif <- function(model_state, X, ...) {
  dif <- model_state$dif
  if (nrow(X) != length(dif$pat))
    stop("A model with `predictors_items` can only be evaluated on the cases ",
         "it was fitted to, since each case's item probabilities depend on its ",
         "own covariate values.", call. = FALSE)
  p     <- pmax(pmin(model_state$parameters$pis, 1 - 1e-15), 1e-15)
  theta <- qlogis(p)
  X0    <- X
  X0[is.na(X0)] <- 0
  M     <- (!is.na(X)) * 1
  out   <- matrix(0, nrow(X), model_state$n_components)
  for (q in seq_len(dif$P)) {
    r   <- dif$rows[[q]]
    eta <- theta + .dif_shift(model_state, q)
    lp  <- plogis(eta, log.p = TRUE)
    lq  <- plogis(-eta, log.p = TRUE)
    out[r, ] <- X0[r, , drop = FALSE] %*% t(lp - lq) +
      M[r, , drop = FALSE] %*% t(lq)
  }
  out
}

#' @exportS3Method
n_parameters.bernoulli_dif <- function(model_state, ...) {
  K    <- model_state$n_components
  free <- model_state$dif$free
  per  <- ifelse(model_state$dif$by_class, K, 1L)
  length(model_state$parameters$pis) + sum(rowSums(free) * per)
}

# The M-step. An item in no named pair gets m_step.bernoulli_nan()'s closed form.
# An item in one gets a weighted binomial regression on the aggregated table:
# one row per (class, covariate pattern) carrying the expected number of
# observed responses and of ones, plus, when the categorical prior is on, one
# row per class at z = 0 carrying the prior's pseudo-counts. Newton from the
# current values with step-halving only ever increases this item's part of the
# expected complete-data objective, so the ECM ascent property holds.
#' @exportS3Method
m_step.bernoulli_dif <- function(model_state, X, resp, weights = NULL,
                                 alpha = NULL, ...) {
  alpha <- alpha %||% .bayes_alpha(model_state, "categorical")
  if (!is.null(weights)) resp <- resp * weights
  dif  <- model_state$dif
  K    <- model_state$n_components
  J    <- ncol(X)
  D    <- ncol(dif$Z)
  prior_obs <- alpha / K
  pis  <- model_state$parameters$pis
  S    <- model_state$parameters$dif_slopes
  if (is.null(S)) S <- matrix(0, K, J * D)

  for (j in seq_len(J)) {
    valid <- !is.na(X[, j])
    if (!any(valid)) next
    marg <- if (is.null(weights)) mean(X[valid, j])
            else sum(X[valid, j] * weights[valid]) / sum(weights[valid])

    if (!any(dif$free[j, ])) {
      rv <- resp[valid, , drop = FALSE]
      pis[, j] <- (colSums(rv * X[valid, j]) + prior_obs * marg) /
        (colSums(rv) + prior_obs)
      next
    }

    n_kp <- s_kp <- matrix(0, K, dif$P)
    for (q in seq_len(dif$P)) {
      r <- dif$rows[[q]]
      r <- r[valid[r]]
      if (!length(r)) next
      rr <- resp[r, , drop = FALSE]
      n_kp[, q] <- colSums(rr)
      s_kp[, q] <- colSums(rr * X[r, j])
    }
    zrows <- dif$Zu
    if (prior_obs > 0) {
      n_kp  <- cbind(n_kp, prior_obs)
      s_kp  <- cbind(s_kp, prior_obs * marg)
      zrows <- rbind(zrows, 0)
    }

    ds  <- which(dif$free[j, ])
    fit <- .dif_item_newton(
      n = as.vector(n_kp), s = as.vector(s_kp),
      A = .dif_item_design(K, zrows, ds, dif$by_class[j]),
      start = c(qlogis(pmax(pmin(pis[, j], 1 - 1e-15), 1e-15)),
                .dif_item_slopes(S, K, J, j, ds, dif$by_class[j])))
    pis[, j] <- plogis(fit[seq_len(K)])
    b <- fit[-seq_len(K)]
    for (i in seq_along(ds)) {
      col <- .dif_col(J, j, ds[i])
      S[, col] <- if (dif$by_class[j]) b[(i - 1L) * K + seq_len(K)] else b[i]
    }
  }
  model_state$parameters$pis        <- pis
  model_state$parameters$dif_slopes <- S
  model_state
}

# Design for one item's regression, rows ordered as as.vector() lays out a
# K x P matrix (class fastest): K threshold columns, then per named covariate
# either one slope column (uniform) or K (one per class).
.dif_item_design <- function(K, zrows, ds, by_class) {
  P  <- nrow(zrows)
  kk <- rep(seq_len(K), P)
  A  <- outer(kk, seq_len(K), "==") * 1
  for (d in ds) {
    z <- rep(zrows[, d], each = K)
    A <- cbind(A, if (by_class) outer(kk, seq_len(K), "==") * z else z)
  }
  A
}

.dif_item_slopes <- function(S, K, J, j, ds, by_class) {
  unlist(lapply(ds, function(d) {
    v <- S[, .dif_col(J, j, d)]
    if (by_class) v else v[1L]
  }), use.names = FALSE)
}

# Newton-Raphson for a binomial logit on aggregated counts, with step-halving so
# the objective never falls. Coefficients are kept within +/-15 on the logit
# scale: past that a probability is within 3e-7 of 0 or 1, so the bound moves
# the log-likelihood by a negligible amount while stopping Newton from chasing
# an estimate that is infinite (a class in which the item is never, or always,
# endorsed).
.dif_item_newton <- function(n, s, A, start, max_iter = 50L, tol = 1e-10) {
  keep <- n > 0
  n <- n[keep]; s <- s[keep]; A <- A[keep, , drop = FALSE]
  obj <- function(b) {
    eta <- drop(A %*% b)
    sum(s * plogis(eta, log.p = TRUE) + (n - s) * plogis(-eta, log.p = TRUE))
  }
  b  <- pmax(pmin(start, 15), -15)
  f0 <- obj(b)
  for (it in seq_len(max_iter)) {
    p <- plogis(drop(A %*% b))
    g <- drop(crossprod(A, s - n * p))
    H <- crossprod(A, A * (n * p * (1 - p)))
    step <- tryCatch(solve(H, g), error = function(e)
      solve(H + diag(1e-8 * max(1, max(diag(H))), ncol(H)), g))
    t <- 1
    repeat {
      b_new <- pmax(pmin(b + t * step, 15), -15)
      f_new <- obj(b_new)
      if (f_new >= f0 - 1e-12 || t < 1e-6) break
      t <- t / 2
    }
    if (f_new < f0 - 1e-12) break
    moved <- max(abs(b_new - b))
    b <- b_new; f0 <- f_new
    if (moved < tol) break
  }
  b
}

# ------------------------------------------------------------------------------
# Packing, for the joint numerical Hessian (R/step3_variance.R)
# ------------------------------------------------------------------------------

.dif_pack <- function(mm) {
  K <- mm$n_components
  J <- ncol(mm$parameters$pis)
  S <- mm$parameters$dif_slopes
  sl <- unlist(lapply(seq_len(J), function(j)
    .dif_item_slopes(S, K, J, j, which(mm$dif$free[j, ]), mm$dif$by_class[j])),
    use.names = FALSE)
  c(qlogis(pmin(pmax(as.vector(mm$parameters$pis), 1e-10), 1 - 1e-10)), sl)
}

.dif_unpack <- function(mm, par) {
  K  <- mm$n_components
  J  <- ncol(mm$parameters$pis)
  mm$parameters$pis[] <- plogis(par[seq_len(K * J)])
  at <- K * J
  for (j in seq_len(J)) {
    for (d in which(mm$dif$free[j, ])) {
      col <- .dif_col(J, j, d)
      if (mm$dif$by_class[j]) {
        mm$parameters$dif_slopes[, col] <- par[at + seq_len(K)]; at <- at + K
      } else {
        mm$parameters$dif_slopes[, col] <- par[at + 1L]; at <- at + 1L
      }
    }
  }
  mm
}

# ------------------------------------------------------------------------------
# The boundary rule
# ------------------------------------------------------------------------------

# A class-specific slope on an item whose probability is pinned at 0 or 1 in a
# class is not identified in that class: the item carries no variation there
# for the covariate to explain. The slope still counts as a parameter, as the
# model defines it, but its value and standard error mean nothing, and saying so
# is the whole remedy -- the same stance R/gaussian_boundary.R takes on a
# collapsed variance.
.check_dif_boundary <- function(fit) {
  mm <- fit$mm
  if (!inherits(mm, "bernoulli_dif")) return(fit)
  items <- colnames(fit$data) %||% paste0("Item", seq_len(ncol(mm$parameters$pis)))
  hits  <- character(0)
  for (j in which(mm$dif$by_class & rowSums(mm$dif$free) > 0)) {
    th <- qlogis(pmin(pmax(mm$parameters$pis[, j], 1e-300), 1 - 1e-16))
    for (k in which(abs(th) >= 9.2))
      hits <- c(hits, sprintf("%s in class %d", items[j], k))
  }
  if (length(hits))
    warning(sprintf(paste0(
      "The class-specific covariate slope of %s is not identified: the item's ",
      "probability in that class is at 0 or 1, so there is no variation for ",
      "the covariate to act on. Its estimate and standard error are ",
      "meaningless; make that item's slope uniform, or drop the pair."),
      paste(hits, collapse = ", ")), call. = FALSE)
  fit
}

# ------------------------------------------------------------------------------
# Reporting
# ------------------------------------------------------------------------------

#' Direct covariate effects on the indicators
#'
#' @description
#' Reports the slopes a model fitted with \code{predictors_items} estimates:
#' the effect of a covariate on the log-odds of endorsing an item, for people
#' in the same latent class. A non-zero slope is differential item
#' functioning (DIF): the item does not measure the classes the same way for
#' everyone. A \emph{uniform} slope is one number shared by every class; a
#' \emph{class-specific} slope has one value per class.
#'
#' The item probabilities the fit reports (\code{measurement_summary()},
#' \code{plot()}) are those of a case whose item covariates are all zero, so
#' centre or code a covariate so that zero is a meaningful reference.
#'
#' @details
#' Which pairs to include is the analyst's choice; the package does not search
#' for them. The usual workflow (Masyn, 2017) fits the model without direct
#' effects, inspects the covariate-item residuals, adds the pairs they point
#' to, and compares the nested fits with \code{\link{lr_test}}. Two cautions
#' from the literature: a DIF screen produces false positives when local
#' independence is also violated, increasingly so at large samples (Depaoli,
#' Jia & Visser, 2025), and a class-specific slope needs considerably larger
#' samples and classes to detect than a uniform one, which is why uniform is
#' the default.
#'
#' Standard errors invert the numerical observed information of the full
#' one-step log-likelihood: item thresholds, direct effects and the
#' class-membership regression together.
#'
#' @param fit A model fitted by \code{fit_mixture()} with \code{predictors_items}.
#' @param se Logical; compute standard errors (costs one numerical Hessian of
#'   the full parameter vector).
#'
#' @return A data frame with one row per slope: \code{item}, \code{covariate},
#'   \code{class} (\code{"all"} for a uniform slope), \code{estimate},
#'   \code{se}, \code{z} and \code{p}.
#'
#' @references
#' Depaoli, S., Jia, F., & Visser, I. (2025). \emph{Structural Equation
#' Modeling}, 32(5), 780-800.
#'
#' Masyn, K. E. (2017). Measurement invariance and differential item
#' functioning in latent class analysis with stepwise multiple indicator
#' multiple cause modeling. \emph{Structural Equation Modeling}, 24(2),
#' 180-197.
#'
#' @seealso \code{\link{fit_mixture}}, \code{\link{lr_test}}
#' @export
dif_effects <- function(fit, se = TRUE) {
  mm <- fit$mm
  if (!inherits(mm, "bernoulli_dif"))
    stop("This model has no `predictors_items`.", call. = FALSE)
  K     <- mm$n_components
  J     <- ncol(mm$parameters$pis)
  items <- colnames(fit$data) %||% paste0("Item", seq_len(J))
  S     <- mm$parameters$dif_slopes

  rows <- list()
  for (j in seq_len(J)) for (d in which(mm$dif$free[j, ])) {
    v   <- S[, .dif_col(J, j, d)]
    cls <- if (mm$dif$by_class[j]) as.character(seq_len(K)) else "all"
    if (!mm$dif$by_class[j]) v <- v[1L]
    rows[[length(rows) + 1L]] <- data.frame(
      item = items[j], covariate = mm$dif$covariates[d], class = cls,
      estimate = v, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$se <- NA_real_

  if (isTRUE(se)) {
    par <- .joint_pack(fit)
    if (!is.null(par)) {
      w  <- fit$sample_weights %||% rep(1, nrow(fit$data))
      ll <- function(v) sum(w * .joint_ll_case(fit, fit$data, v))
      V  <- pinv(-.step1_fd_hessian(ll, par))
      n_sl <- nrow(out)
      at   <- K * J
      out$se <- sqrt(pmax(diag(V)[at + seq_len(n_sl)], 0))
    }
  }
  out$z <- out$estimate / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  rownames(out) <- NULL
  out
}
