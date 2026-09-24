# ==============================================================================
# A distal outcome on the last occasion of the three-step latent transition
# model.
# ==============================================================================
#
# Step 3 of the three-step (R/lta_threestep.R) models the assigned statuses
# W_1..W_T as indicators of the true ones with a known classification error.
# A distal outcome measured after the last occasion adds one factor to that
# occasion's emission,
#
#   P(W_T | S_T) * prod_p f(Y_p | S_T),
#
# where each Y_p is either Gaussian with a status-specific mean and variance or
# binary with a status-specific probability. The factor is estimated jointly
# with the structural model, the ML correction of Vermunt (2010) with the
# outcome as one more term of step 3's likelihood: it moves the posterior of S_T and
# with it the transitions, so fitting the two one after the other would give a
# different answer.
#
# The outcome is carried as extra trailing columns of the data matrix the
# recursion reads, and its parameters on the measurement model as `mm$distal`.
# That puts it everywhere an emission is formed -- the E-step, the case
# log-likelihood behind the Hessian, the posteriors a fit reports -- and keeps
# it row-aligned with the data through every subset the variance code takes,
# without a second argument threaded through each of them. The assigned-status
# emission stays frozen; only the outcome's parameters move in the M-step.
#
# Modal assignment only. Under proportional assignment step 3's rows are
# status combinations with no case of their own, so there is nothing for a
# case's outcome to attach to.

# The user's outcome -> a numeric matrix plus a type per column. A column is
# binary when its observed values are 0/1 (or logical, or a two-level factor,
# coded 1 for the second level); anything else numeric is Gaussian.
.lta_distal_prepare <- function(distal, n) {
  if (is.null(distal)) return(NULL)
  if (is.vector(distal) || is.factor(distal)) distal <- data.frame(distal)
  distal <- as.data.frame(distal)
  if (nrow(distal) != n)
    stop(sprintf(paste0(
      "`distal` has %d rows and the indicators have %d. It must hold one row ",
      "per case, in the same order."), nrow(distal), n), call. = FALSE)
  if (ncol(distal) < 1L)
    stop("`distal` has no columns.", call. = FALSE)
  nm <- names(distal)
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("distal", seq_len(ncol(distal)))

  type <- character(ncol(distal))
  Y <- matrix(NA_real_, n, ncol(distal), dimnames = list(NULL, nm))
  for (p in seq_len(ncol(distal))) {
    v <- distal[[p]]
    if (is.factor(v) || is.character(v)) {
      v <- factor(v)
      if (nlevels(v) != 2L)
        stop(sprintf(paste0(
          "Distal outcome `%s` is categorical with %d categories; only ",
          "binary and continuous outcomes are supported."), nm[p], nlevels(v)),
          call. = FALSE)
      Y[, p] <- as.integer(v) - 1L
      type[p] <- "binary"
    } else if (is.logical(v) || is.numeric(v)) {
      v <- as.numeric(v)
      obs <- v[!is.na(v)]
      if (length(obs) < 2L)
        stop(sprintf("Distal outcome `%s` has fewer than two observed values.",
                     nm[p]), call. = FALSE)
      Y[, p] <- v
      type[p] <- if (all(obs %in% c(0, 1))) "binary" else "gaussian"
    } else {
      stop(sprintf("Distal outcome `%s` must be numeric, logical or a factor.",
                   nm[p]), call. = FALSE)
    }
  }
  list(Y = Y, type = type, names = nm)
}

# Starting values: the outcome's moments within each assigned status at the
# last occasion. The classification error is small wherever the three-step is
# worth running, so these sit close to the answer, and EM does the rest.
.lta_distal_init <- function(d, assigned, w, K) {
  P  <- length(d$type)
  mu <- v <- pr <- matrix(NA_real_, K, P, dimnames = list(NULL, d$names))
  for (p in seq_len(P)) {
    y  <- d$Y[, p]
    ok <- !is.na(y)
    tot_m <- stats::weighted.mean(y[ok], w[ok])
    tot_v <- sum(w[ok] * (y[ok] - tot_m)^2) / sum(w[ok])
    for (k in seq_len(K)) {
      i <- ok & assigned == k
      m <- if (any(i)) stats::weighted.mean(y[i], w[i]) else tot_m
      if (d$type[p] == "gaussian") {
        mu[k, p] <- m
        s2 <- if (sum(i) > 1L) sum(w[i] * (y[i] - m)^2) / sum(w[i]) else tot_v
        v[k, p]  <- max(s2, 1e-3 * tot_v)
      } else {
        pr[k, p] <- min(max(m, 0.01), 0.99)
      }
    }
  }
  list(type = d$type, names = d$names, mean = mu, var = v, prob = pr)
}

# The outcome's columns in the data matrix: everything after the indicators.
.lta_distal_cols <- function(mm) mm$n_items * mm$n_times + seq_along(mm$distal$type)

# n x K log-density of the outcome given each status. A missing value
# contributes nothing, which is FIML on the outcome.
.lta_distal_loglik <- function(dist, Y) {
  K   <- nrow(if (is.null(dist$mean)) dist$prob else dist$mean)
  out <- matrix(0, nrow(Y), K)
  for (p in seq_along(dist$type)) {
    y  <- Y[, p]
    ok <- !is.na(y)
    for (k in seq_len(K)) {
      out[ok, k] <- out[ok, k] + if (dist$type[p] == "gaussian")
        stats::dnorm(y[ok], dist$mean[k, p], sqrt(dist$var[k, p]), log = TRUE)
      else
        stats::dbinom(y[ok], 1L, dist$prob[k, p], log = TRUE)
    }
  }
  out
}

# Weighted maximum likelihood given the last occasion's posteriors. No prior:
# a status-specific variance can in principle collapse onto one case, so it is
# floored at a millionth of the outcome's total variance, far below anything a
# real solution reports.
.lta_distal_mstep <- function(dist, Y, g, w) {
  for (p in seq_along(dist$type)) {
    y  <- Y[, p]
    ok <- !is.na(y)
    R  <- g[ok, , drop = FALSE] * w[ok]
    nk <- colSums(R)
    m  <- as.vector(crossprod(R, y[ok])) / pmax(nk, 1e-300)
    if (dist$type[p] == "gaussian") {
      tot <- sum(w[ok] * (y[ok] - stats::weighted.mean(y[ok], w[ok]))^2) /
        sum(w[ok])
      dist$mean[, p] <- m
      dist$var[, p]  <- pmax(colSums(R * outer(y[ok], m, "-")^2) /
                               pmax(nk, 1e-300), 1e-6 * tot)
    } else {
      dist$prob[, p] <- pmin(pmax(m, 1e-10), 1 - 1e-10)
    }
  }
  dist
}

.lta_distal_n_params <- function(dist) {
  if (is.null(dist)) return(0L)
  K <- nrow(if (is.null(dist$mean)) dist$prob else dist$mean)
  as.integer(K * sum(ifelse(dist$type == "gaussian", 2L, 1L)))
}

# The outcome's blocks in .lta_par_layout()'s vector: per outcome, the K means
# then the K log standard deviations, or the K logits of P(Y = 1).
.lta_distal_layout <- function(dist, K) {
  out <- list()
  for (p in seq_along(dist$type)) {
    if (dist$type[p] == "gaussian") {
      out[[length(out) + 1L]] <- list(kind = "distal_mu", p = p, len = K)
      out[[length(out) + 1L]] <- list(kind = "distal_log_sd", p = p, len = K)
    } else {
      out[[length(out) + 1L]] <- list(kind = "distal_logit", p = p, len = K)
    }
  }
  out
}

# The per-status table a user reads, from the fitted outcome and the rows of
# the covariance its blocks occupy. Standard errors on the variance and the
# probability are carried from the packed scale by the delta method.
.lta_distal_table <- function(fit) {
  dist <- fit$mm$distal
  K    <- fit$n_statuses
  lay  <- .lta_par_layout(fit)
  len  <- vapply(lay, function(b) as.integer(b$len), integer(1))
  beg  <- cumsum(len) - len
  V    <- fit$se$vcov
  se_of <- function(kind, p) {
    i <- which(vapply(lay, function(b)
      identical(b$kind, kind) && identical(b$p, p), logical(1)))
    if (is.null(V) || length(i) != 1L || nrow(V) < beg[i] + K)
      return(rep(NA_real_, K))
    sqrt(pmax(diag(V)[beg[i] + seq_len(K)], 0))
  }
  rows <- list()
  for (p in seq_along(dist$type)) {
    if (dist$type[p] == "gaussian") {
      rows[[length(rows) + 1L]] <- data.frame(
        outcome = dist$names[p], status = seq_len(K), parameter = "mean",
        estimate = dist$mean[, p], se = se_of("distal_mu", p))
      rows[[length(rows) + 1L]] <- data.frame(
        outcome = dist$names[p], status = seq_len(K), parameter = "variance",
        estimate = dist$var[, p],
        se = 2 * dist$var[, p] * se_of("distal_log_sd", p))
    } else {
      pr <- dist$prob[, p]
      rows[[length(rows) + 1L]] <- data.frame(
        outcome = dist$names[p], status = seq_len(K), parameter = "P(Y = 1)",
        estimate = pr, se = pr * (1 - pr) * se_of("distal_logit", p))
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Status-vs-status contrasts on one outcome, for outcome_contrasts(): the
# difference in means for a Gaussian outcome, in log odds of Y = 1 for a
# binary one, with the covariance the fit carries.
.lta_distal_contrast_parts <- function(fit, outcome = NULL) {
  dist <- fit$mm$distal
  if (is.null(dist))
    stop("This latent transition model has no distal outcome. Fit one with ",
         "`fit_lta(..., n_steps = 3, distal = ...)`.", call. = FALSE)
  p <- if (is.null(outcome)) 1L else if (is.character(outcome))
    match(outcome, dist$names) else as.integer(outcome)
  if (length(p) != 1L || is.na(p) || p < 1L || p > length(dist$type))
    stop("`outcome` must name one of the distal outcomes: ",
         paste(dist$names, collapse = ", "), ".", call. = FALSE)
  if (is.null(fit$se$vcov))
    stop("The fit carries no covariance matrix; refit with ",
         "`standard_errors = TRUE`.", call. = FALSE)
  K    <- fit$n_statuses
  lay  <- .lta_par_layout(fit)
  len  <- vapply(lay, function(b) as.integer(b$len), integer(1))
  beg  <- cumsum(len) - len
  kind <- if (dist$type[p] == "gaussian") "distal_mu" else "distal_logit"
  i    <- which(vapply(lay, function(b)
    identical(b$kind, kind) && identical(b$p, p), logical(1)))
  cols <- beg[i] + seq_len(K)
  theta <- if (dist$type[p] == "gaussian") dist$mean[, p] else
    stats::qlogis(dist$prob[, p])
  list(type   = if (dist$type[p] == "gaussian") "continuous" else "categorical",
       method = fit$se$method %||% "observed information",
       label  = dist$names[p],
       parts  = list(list(category = if (dist$type[p] == "gaussian")
                            NA_integer_ else 1L,
                          theta = theta,
                          V = fit$se$vcov[cols, cols, drop = FALSE],
                          index = seq_len(K))))
}
