# The reproducibility contract of R/parallel.R: a fit is the same at every
# `n_cores`, because the workers never draw random numbers. Every number the
# validation suite asserts now comes off a parallel run, so that equality is
# what makes the suite's targets mean anything.
#
# Every test here spawns worker processes, so all of them skip on CRAN. The
# worker count is 2 rather than more for the same reason R CMD check caps it
# there (see .par_max_cores()).

.par_sim <- function(seed = 11, n = 300, K = 3, J = 6) {
  set.seed(seed)
  pis <- matrix(runif(K * J, 0.15, 0.85), K, J)
  cls <- sample(K, n, TRUE)
  matrix(rbinom(n * J, 1, pis[cls, ]), n, J)
}

.par_fit <- function(X, ...) suppressMessages(suppressWarnings(
  fit_mixture(X, n_classes = 3, measurement = "binary",
              n_init = 6, random_state = 7, ...)))

test_that(".default_n_cores() reads the option and refuses nonsense", {
  old <- getOption("mixtureEM.n_cores")
  on.exit(options(mixtureEM.n_cores = old), add = TRUE)

  options(mixtureEM.n_cores = NULL)
  expect_identical(.default_n_cores(), 1L)

  options(mixtureEM.n_cores = 4)
  expect_identical(.default_n_cores(), 4L)

  # A bad value falls back to one worker rather than erroring: this is read on
  # every fit, and a typo in an option set weeks ago should not stop a model
  # from being estimated.
  for (bad in list("junk", 0, -3, NA, c(2, 3))) {
    options(mixtureEM.n_cores = bad)
    expect_identical(.default_n_cores(), 1L)
  }
})

test_that("a fit is identical sequentially and on workers", {
  skip_on_cran()
  on.exit(.par_stop_cluster(), add = TRUE)

  X   <- .par_sim()
  seq <- .par_fit(X, n_cores = 1L)
  par <- .par_fit(X, n_cores = 2L)

  # Not `tolerance =`: the claim is bit-for-bit, not "close enough".
  expect_identical(par$metrics$ll, seq$metrics$ll)
  expect_identical(par$start_lls,  seq$start_lls)
  expect_identical(par$log_resp,   seq$log_resp)
  expect_identical(par$weights,    seq$weights)
  expect_identical(par$mm$parameters$pis, seq$mm$parameters$pis)
})

test_that("the option supplies the default and an argument overrides it", {
  skip_on_cran()
  old <- getOption("mixtureEM.n_cores")
  on.exit({ options(mixtureEM.n_cores = old); .par_stop_cluster() }, add = TRUE)

  X   <- .par_sim()
  seq <- .par_fit(X, n_cores = 1L)

  options(mixtureEM.n_cores = 2)
  expect_identical(.par_fit(X)$metrics$ll, seq$metrics$ll)
  expect_identical(.par_state$n, 2L)

  # An explicit argument wins over the option, and must not disturb the pool.
  .par_stop_cluster()
  expect_identical(.par_fit(X, n_cores = 1L)$metrics$ll, seq$metrics$ll)
  expect_null(.par_state$cl)
})

test_that("the worker pool is built once and reused", {
  skip_on_cran()
  on.exit(.par_stop_cluster(), add = TRUE)
  .par_stop_cluster()

  X <- .par_sim()
  .par_fit(X, n_cores = 2L)
  first <- .par_state$cl
  expect_false(is.null(first))

  # The second fit must not rebuild: paying the spawn and the package load per
  # call is what made `n_cores` a pessimisation for anything but the largest
  # jobs, and is the whole reason the pool is cached.
  .par_fit(X, n_cores = 2L)
  expect_identical(.par_state$cl, first)

  # A job with fewer pieces than workers reuses the pool rather than resizing
  # it. Sizing on length(X) would rebuild whenever `n_init` changed.
  suppressMessages(suppressWarnings(
    fit_mixture(X, n_classes = 3, measurement = "binary",
                n_init = 2, random_state = 7, n_cores = 2L)))
  expect_identical(.par_state$cl, first)
})

test_that(".par_lapply() falls back to sequential rather than failing", {
  expect_identical(.par_lapply(1:3, function(i) i * 2, n_cores = 1L),
                   list(2, 4, 6))
  # One element: not worth a worker, and the early return must still be a list.
  expect_identical(.par_lapply(list(5), function(i) i + 1, n_cores = 4L),
                   list(6))
})

test_that("a dead pool is detected and rebuilt", {
  skip_on_cran()
  on.exit(.par_stop_cluster(), add = TRUE)

  expect_false(is.null(.par_cluster(2L)))
  dead <- .par_state$cl

  # Kill the workers behind the cached handle without clearing the cache, which
  # is what a crashed worker looks like from here.
  suppressWarnings(try(parallel::stopCluster(dead), silent = TRUE))
  expect_false(.par_cluster_alive(dead))

  rebuilt <- .par_cluster(2L)
  expect_false(is.null(rebuilt))
  expect_false(identical(rebuilt, dead))

  # The worker gets a function whose environment is the package namespace when
  # it is written inline here, and serialising that warns. Package code never
  # hands over such a closure; a test should not either.
  square <- function(i) i^2
  environment(square) <- baseenv()
  expect_identical(unlist(.par_lapply(1:4, square, n_cores = 2L)),
                   c(1, 4, 9, 16))
})
