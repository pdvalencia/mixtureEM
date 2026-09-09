test_that("the scaled recursion returns exactly what the log-space one returns", {
  # The scaled forward-backward is the same arithmetic in the probability
  # domain, so this is an equality check, not a tolerance-shopping exercise: any
  # disagreement above floating-point noise is a defect in one of them. The
  # log-space form is the reference because it is the one every locked target in
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

test_that("a covariate-transition model still takes the log-space recursion", {
  # The dispatcher's whole job. If this ever silently routed a list-valued
  # `log_tau` into the scaled form, the matrix product would be wrong rather
  # than an error, so the branch is asserted rather than assumed.
  set.seed(7)
  n <- 60; K <- 3L; Tn <- 2L; J <- 4L
  logB <- lapply(seq_len(Tn), function(t) {
    p <- matrix(runif(K * J, 0.1, 0.9), K, J)
    X <- matrix(rbinom(n * J, 1L, 0.5), n, J)
    X %*% t(log(p)) + (1 - X) %*% t(log1p(-p))
  })
  d <- runif(K); log_delta <- log(d / sum(d))
  log_tau <- list(lapply(seq_len(K), function(k) {
    m <- matrix(runif(n * K), n, K); log(m / rowSums(m))
  }))
  w <- rep(1, n)

  got <- .lta_forward_backward(logB, log_delta, log_tau, w, keep_pairwise = TRUE)
  ref <- .lta_fb_log(logB, log_delta, log_tau, w, keep_pairwise = TRUE)
  expect_equal(got$ll, ref$ll)
  expect_equal(got$xi[[1]], ref$xi[[1]])
})
