test_that("the scaled recursion returns exactly what the log-space one returns", {
  # The scaled forward-backward is the same arithmetic in the probability
  # domain, so this is an equality check, not a tolerance-shopping exercise: any
  # disagreement above floating-point noise is a defect in one of them. The
  # log-space form is the oracle because it is the one every locked benchmark in
  # the package was measured on.
  set.seed(202)

  check <- function(n, K, Tn, J, cov_delta = FALSE, extreme = FALSE) {
    logB <- lapply(seq_len(Tn), function(t) {
      p <- matrix(runif(K * J, if (extreme) 0.001 else 0.05,
                        if (extreme) 0.999 else 0.95), K, J)
      X <- matrix(rbinom(n * J, 1L, 0.5), n, J)
      X %*% t(log(p)) + (1 - X) %*% t(log1p(-p))
    })
    log_delta <- if (cov_delta) {
      d <- matrix(runif(n * K), n, K); log(d / rowSums(d))
    } else {
      d <- runif(K); log(d / sum(d))
    }
    log_tau <- lapply(seq_len(max(Tn - 1L, 0L)), function(t) {
      m <- matrix(runif(K * K), K, K)
      log(m / rowSums(m))
    })
    w <- runif(n, 0.5, 2)

    a <- .lta_fb_scaled(logB, log_delta, log_tau, w, keep_pairwise = TRUE)
    b <- .lta_fb_log(logB, log_delta, log_tau, w, keep_pairwise = TRUE)

    expect_equal(a$ll, b$ll, tolerance = 1e-10)
    for (t in seq_len(Tn)) expect_equal(a$gamma[[t]], b$gamma[[t]],
                                        tolerance = 1e-10)
    for (t in seq_len(max(Tn - 1L, 0L))) {
      expect_equal(a$xi[[t]], b$xi[[t]], tolerance = 1e-10)
      for (k in seq_len(K))
        expect_equal(a$pairwise[[t]][[k]], b$pairwise[[t]][[k]],
                     tolerance = 1e-10)
    }
  }

  check(n = 200, K = 3, Tn = 3, J = 5)
  check(n = 150, K = 4, Tn = 2, J = 6)
  check(n = 120, K = 2, Tn = 5, J = 4)
  check(n = 100, K = 3, Tn = 1, J = 5)            # no transitions at all
  check(n = 100, K = 3, Tn = 3, J = 5, cov_delta = TRUE)

  # Forty items with response probabilities out at 0.001 puts the log emission
  # near -280 for some cases, which is where a probability-domain recursion is
  # usually said to fail. It does not here, because each occasion's row maximum
  # is factored out before anything is exponentiated.
  check(n = 200, K = 3, Tn = 3, J = 40, extreme = TRUE)
})

test_that("case-varying transitions give the log-space answer through the scaled recursion", {
  # A covariate-transition model hands a list over origins of n x K matrices.
  # The scaled form writes each matrix product row by row for it; the log-space
  # form stays the oracle, and is what answers when a case underflows.
  set.seed(7)
  check <- function(n, K, Tn, J, cov_delta = FALSE, zero = FALSE) {
    logB <- lapply(seq_len(Tn), function(t) {
      p <- matrix(runif(K * J, 0.1, 0.9), K, J)
      X <- matrix(rbinom(n * J, 1L, 0.5), n, J)
      X %*% t(log(p)) + (1 - X) %*% t(log1p(-p))
    })
    log_delta <- if (cov_delta) {
      d <- matrix(runif(n * K), n, K); log(d / rowSums(d))
    } else {
      d <- runif(K); log(d / sum(d))
    }
    log_tau <- lapply(seq_len(Tn - 1L), function(t)
      lapply(seq_len(K), function(k) {
        m <- matrix(runif(n * K), n, K)
        if (zero) m[, (k %% K) + 1L] <- 0           # a forbidden transition
        log(m / rowSums(m))
      }))
    w <- runif(n, 0.5, 2)
    got <- .lta_forward_backward(logB, log_delta, log_tau, w, keep_pairwise = TRUE)
    ref <- .lta_fb_log(logB, log_delta, log_tau, w, keep_pairwise = TRUE)
    expect_equal(got$ll, ref$ll, tolerance = 1e-10)
    for (t in seq_len(Tn)) expect_equal(got$gamma[[t]], ref$gamma[[t]],
                                        tolerance = 1e-10)
    for (t in seq_len(Tn - 1L)) {
      expect_equal(got$xi[[t]], ref$xi[[t]], tolerance = 1e-10)
      for (k in seq_len(K))
        expect_equal(got$pairwise[[t]][[k]], ref$pairwise[[t]][[k]],
                     tolerance = 1e-10)
    }
  }
  check(n = 60, K = 3, Tn = 2, J = 4)
  check(n = 80, K = 4, Tn = 3, J = 5, cov_delta = TRUE)
  check(n = 50, K = 3, Tn = 4, J = 4, zero = TRUE)
  check(n = 1, K = 3, Tn = 2, J = 4)

  # One case whose every path runs through a probability of exp(-800): the
  # scaled alpha underflows to zero, the log-space one does not, and the pass
  # falls back to it.
  logB    <- list(matrix(0, 1, 2), matrix(c(-800, 0), 1, 2))
  log_tau <- list(list(matrix(c(0, -800), 1, 2), matrix(log(0.5), 1, 2)))
  expect_false(is.finite(.lta_fb_scaled(logB, c(0, -800), log_tau, 1)$ll))
  got <- .lta_forward_backward(logB, c(0, -800), log_tau, 1, keep_pairwise = TRUE)
  expect_equal(got$ll, -800 + log(2.5))
  expect_equal(got$gamma[[2]], matrix(c(0.4, 0.6), 1, 2))
})
