# ==============================================================================
# Latent transition analysis - global (Viterbi) decoding
# ==============================================================================
#
# The oracle here is exhaustive enumeration of all K^T status sequences, so
# these tests need no reference program: for the small K and T used below the
# brute-force maximum IS the answer the recursion must reproduce.

# Enumerate every status sequence and return the per-case maximiser.
.brute_force_path <- function(logB, log_delta, log_tau) {
  Tn <- length(logB); n <- nrow(logB[[1]]); K <- ncol(logB[[1]])
  grid <- as.matrix(expand.grid(rep(list(seq_len(K)), Tn)))
  best <- matrix(0L, n, Tn); bestv <- rep(-Inf, n)
  for (r in seq_len(nrow(grid))) {
    s <- grid[r, ]
    v <- if (is.matrix(log_delta)) log_delta[, s[1]] else rep(log_delta[s[1]], n)
    v <- v + logB[[1]][, s[1]]
    if (Tn > 1L) for (t in 2:Tn) {
      LT <- log_tau[[t - 1]]
      tr <- if (is.list(LT)) LT[[s[t - 1]]][, s[t]] else LT[s[t - 1], s[t]]
      v <- v + tr + logB[[t]][, s[t]]
    }
    hit <- v > bestv
    if (any(hit)) {
      bestv[hit] <- v[hit]
      best[hit, ] <- matrix(s, sum(hit), Tn, byrow = TRUE)
    }
  }
  list(path = best, logp = bestv)
}

.fit_pieces <- function(fit, class = NULL) {
  st <- if (is.null(class)) fit else mixtureEM:::.lta_class_state(fit, class)
  list(logB = mixtureEM:::.lta_emission_loglik(fit$mm, fit$data),
       log_delta = mixtureEM:::.lta_log_delta(st),
       log_tau   = mixtureEM:::.lta_log_tau(st))
}


# Strip every attribute except dim, so a path matrix can be compared with a
# bare brute-force one: unname() leaves "probability" and "class_assigned".
.bare <- function(m) { attributes(m) <- list(dim = dim(m)); m }

.sim_lta <- function(seed, n, Tn, J, K, sep, stay) {
  set.seed(seed)
  S <- matrix(0L, n, Tn); S[, 1] <- sample(seq_len(K), n, TRUE)
  for (t in 2:Tn)
    S[, t] <- ifelse(stats::runif(n) < stay, S[, t - 1], sample(seq_len(K), n, TRUE))
  p <- matrix(0, K, J)
  for (k in seq_len(K)) p[k, ] <- (1 - sep) / 2 + sep * (k - 1) / max(K - 1, 1)
  X <- matrix(0, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J))
    X[, (t - 1) * J + j] <- stats::rbinom(n, 1, p[S[, t], j])
  X
}

# ------------------------------------------------------------------------------
# The recursion is exactly the enumeration
# ------------------------------------------------------------------------------

test_that("viterbi decoding reproduces exhaustive enumeration", {
  X <- .sim_lta(31, 400, 4, 6, 3, 0.94, 0.8)
  fit <- fit_lta(X, n_statuses = 3, times = 4, measurement = "binary",
                 n_init = 10, random_state = 31, standard_errors = FALSE)
  p <- .fit_pieces(fit)
  b <- .brute_force_path(p$logB, p$log_delta, p$log_tau)
  v <- class_assignments(fit, "viterbi")

  expect_equal(.bare(v), b$path)
  expect_equal(dim(v), c(400L, 4L))
  expect_equal(colnames(v), fit$longitudinal$time_labels)
  expect_equal(unname(attr(v, "probability")), exp(b$logp - fit$ll_case),
               tolerance = 1e-10)
  # A path probability is a probability.
  expect_true(all(attr(v, "probability") > 0 & attr(v, "probability") <= 1))
})

test_that("viterbi decoding matches enumeration under a covariate on the initial status", {
  X <- .sim_lta(34, 500, 4, 4, 2, 0.8, 0.7)
  z <- stats::rnorm(500)
  fit <- fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
                 predictors_initial = z, n_init = 8, random_state = 35,
                 standard_errors = FALSE)
  p <- .fit_pieces(fit)
  # A covariate on the initial status makes log_delta an n x K matrix.
  expect_true(is.matrix(p$log_delta))
  b <- .brute_force_path(p$logB, p$log_delta, p$log_tau)
  expect_equal(.bare(class_assignments(fit, "viterbi")), b$path)
})

test_that("viterbi decoding matches enumeration in a multiple-group model", {
  X <- .sim_lta(34, 500, 4, 4, 2, 0.8, 0.7)
  g <- rep(c("a", "b"), each = 250)
  fit <- fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
                 group = g, group_effects = "both", n_init = 8,
                 random_state = 34, standard_errors = FALSE)
  p <- .fit_pieces(fit)
  # A group enters as a covariate, so the transitions are per case, not per
  # group: each element of log_tau is a list of K matrices of size n x K.
  expect_true(is.list(p$log_tau[[1]]))
  expect_equal(dim(p$log_tau[[1]][[1]]), c(500L, 2L))
  b <- .brute_force_path(p$logB, p$log_delta, p$log_tau)
  expect_equal(.bare(class_assignments(fit, "viterbi")), b$path)
})

test_that("a fully missing occasion contributes only the transition term", {
  X <- .sim_lta(101, 200, 4, 3, 2, 0.7, 0.75)
  X[1:25, mixtureEM:::.time_block_cols(3, 3)] <- NA
  fit <- fit_lta(X, n_statuses = 2, times = 4, measurement = "binary",
                 n_init = 5, random_state = 8, standard_errors = FALSE)
  p <- .fit_pieces(fit)
  expect_true(all(p$logB[[3]][1:25, ] == 0))
  b <- .brute_force_path(p$logB, p$log_delta, p$log_tau)
  expect_equal(.bare(class_assignments(fit, "viterbi")), b$path)
})

# ------------------------------------------------------------------------------
# The reason the feature exists
# ------------------------------------------------------------------------------

test_that("viterbi never steps through a forbidden transition", {
  K <- 3; Tn <- 5; n <- 400
  set.seed(22)
  S <- matrix(0L, n, Tn); S[, 1] <- sample(1:2, n, TRUE, c(0.7, 0.3))
  for (t in 2:Tn) S[, t] <- pmin(K, S[, t - 1] + stats::rbinom(n, 1, 0.4))
  p <- matrix(0, K, 2)
  for (k in 1:K) p[k, ] <- 0.30 + (k - 1) * 0.40 / (K - 1)
  X <- matrix(0, n, Tn * 2)
  for (t in seq_len(Tn)) for (j in 1:2)
    X[, (t - 1) * 2 + j] <- stats::rbinom(n, 1, p[S[, t], j])
  forb <- matrix(0, K, K); forb[lower.tri(forb)] <- 1   # no returning to an earlier stage

  fit <- suppressWarnings(
    fit_lta(X, n_statuses = K, times = Tn, measurement = "binary",
            forbidden_transitions = forb, n_init = 10, random_state = 22,
            standard_errors = FALSE))

  crossings <- function(path) {
    out <- rep(FALSE, nrow(path))
    for (t in seq_len(ncol(path) - 1L))
      out <- out | (fit$tau[[t]][cbind(path[, t], path[, t + 1L])] == 0)
    out
  }
  # The guarantee. A path the model gives probability zero can never be the
  # most probable path, so global decoding cannot return one.
  expect_equal(sum(crossings(class_assignments(fit, "viterbi"))), 0L)
  # Local decoding carries no such guarantee; this is what the feature is for.
  # Not asserted as a positive count - it depends on how well separated the
  # statuses are - but the two decodings must genuinely differ here.
  expect_gt(sum(rowSums(class_assignments(fit) !=
                          class_assignments(fit, "viterbi")) > 0), 0L)
})

# ------------------------------------------------------------------------------
# A mixture over chains
# ------------------------------------------------------------------------------

test_that("a mixture is decoded jointly over class and path", {
  set.seed(61); n <- 500; Tn <- 4; J <- 3; noise <- 0.32
  stayer <- stats::rbinom(n, 1, .5) == 1
  S <- matrix(0L, n, Tn); S[, 1] <- stats::rbinom(n, 1, .5) + 1L
  for (t in 2:Tn)
    S[, t] <- ifelse(stayer, S[, t - 1],
                     ifelse(stats::runif(n) < .5, S[, t - 1], 3L - S[, t - 1]))
  pr <- rbind(rep(noise, J), rep(1 - noise, J))
  X <- matrix(0, n, Tn * J)
  for (t in seq_len(Tn)) for (j in seq_len(J))
    X[, (t - 1) * J + j] <- stats::rbinom(n, 1, pr[S[, t], j])
  fit <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = Tn, measurement = "binary",
            mover_stayer = TRUE, n_init = 15, random_state = 61,
            standard_errors = FALSE))
  skip_if_not(identical(fit$n_classes, 2L))

  v   <- class_assignments(fit, "viterbi")
  cls <- attr(v, "class_assigned")
  expect_equal(length(cls), n)
  expect_true(all(cls %in% 1:2))

  # Enumerate over (class, path) jointly.
  best <- matrix(0L, n, Tn); bestv <- rep(-Inf, n); bestc <- integer(n)
  for (c in 1:2) {
    p  <- .fit_pieces(fit, c)
    b  <- .brute_force_path(p$logB, p$log_delta, p$log_tau)
    lv <- b$logp + log(fit$class_weights[c])
    h  <- lv > bestv
    bestv[h] <- lv[h]; best[h, ] <- b$path[h, ]; bestc[h] <- c
  }
  expect_equal(.bare(v), best)
  expect_equal(cls, bestc)

  # The stayer is the last class and its only admissible move is to stay put,
  # so a case decoded into it must have a constant path. Local decoding, which
  # reads the class-mixed posterior one occasion at a time, does not know that.
  stayers <- cls == 2L
  expect_equal(sum(stayers & rowSums(v != v[, 1]) > 0), 0L)
})
