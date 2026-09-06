# ==============================================================================
# Optional parallelism for the repeated-fit loops
# ==============================================================================
#
# Several operations in the package are a loop of independent model fits:
# multi-start estimation, the bootstrap likelihood ratio test, the bootstrap
# standard errors, and the K grids of the compare_*() functions. A BLRT at the
# defaults is over two thousand fits. Nothing in those loop bodies depends on
# any other iteration, so they are worth spreading over cores.
#
# `parallel` is a base R package, so this adds no third-party dependency.
#
# REPRODUCIBILITY
#
# Workers never draw random numbers. Every random quantity a loop needs -- the
# starting values of a restart, the synthetic data of a bootstrap replicate --
# is drawn in the parent process, in the same order the sequential code drew it,
# and handed to the worker as data. What the worker does with it is a
# deterministic function of its input.
#
# The consequence is the property the validation suite depends on: for a given
# seed the results are identical for every value of `n_cores`, including the
# sequential path, so turning cores on can never move a published number.

# Largest worker count we will actually start. R CMD check sets
# _R_CHECK_LIMIT_CORES_, and CRAN's policy caps checks at two cores.
.par_max_cores <- function() {
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && !identical(tolower(chk), "false")) return(2L)
  n <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
  if (is.na(n)) 1L else max(1L, n)
}

# Path to this package's source tree when it was loaded by devtools/pkgload
# rather than installed, otherwise NULL. An installed package carries a Meta/
# directory; a source tree under load_all() does not. Workers are fresh R
# processes, so under load_all() they would otherwise pick up whatever release
# happens to be installed instead of the code under test.
.par_dev_path <- function() {
  p <- tryCatch(getNamespaceInfo(asNamespace("mixtureEM"), "path"),
                error = function(e) NULL)
  if (is.null(p) || !nzchar(p)) return(NULL)
  if (dir.exists(file.path(p, "Meta"))) return(NULL)
  if (!file.exists(file.path(p, "DESCRIPTION"))) return(NULL)
  p
}

# The session-wide default for `n_cores`, read by every user-facing fit.
#
# The shipped default is one worker. Parallelism costs a PSOCK cluster and a
# copy of the data per worker, which is a net loss on a fit that takes a
# second, and a package must not start processes uninvited during a check.
#
# `options(mixtureEM.n_cores = )` moves it for a whole session, which is what a
# long analysis wants: a value passed to one call is forgotten by the next, and
# the repeated-fit functions -- blrt(), compare_mixtures(),
# bootstrap_covariates(), bivariate_residuals() -- are precisely the slow ones
# a user forgets to pass it to. An argument given explicitly always wins over
# the option.
#
# Nothing about reproducibility changes. Workers never draw random numbers (see
# the note at the top of this file), so a fit is identical at any worker count,
# and the option cannot move a published number.
.default_n_cores <- function() {
  n <- suppressWarnings(as.integer(getOption("mixtureEM.n_cores", 1L)))
  if (length(n) != 1L || is.na(n) || n < 1L) 1L else n
}

# ------------------------------------------------------------------------------
# The worker pool, built once and kept
# ------------------------------------------------------------------------------
#
# A cluster used to be built and torn down inside every .par_lapply() call. That
# is affordable when the parallel call is the whole job -- a BLRT's two thousand
# fits -- and ruinous when it is one fit among many, which is what a validation
# run, a compare_mixtures() sweep or an interactive session actually looks like.
#
# Measured on this machine, under pkgload::load_all(), building the pool and
# loading the package onto it:
#
#    4 workers    spawn 0.18s   load package 1.93s   TOTAL 2.12s
#   14 workers    spawn 0.39s   load package 10.44s  TOTAL 10.86s
#
# Ten seconds per call is more than most single fits cost, so paying it per call
# made `n_cores` a pessimisation for every job except the very largest. The pool
# is therefore cached for the session and reused. An installed package loads far
# faster than a load_all() source tree, so the figures above are the worst case,
# but the argument holds either way.
#
# Nothing here is created unless a caller asks for more than one worker, so a
# default session -- and R CMD check -- never starts a process at all.
.par_state <- new.env(parent = emptyenv())

# Is this cluster still usable? A worker can die between calls, and a stale
# handle fails inside parLapply() where the fit is already underway. One
# round-trip is cheap beside what it protects.
.par_cluster_alive <- function(cl) {
  if (is.null(cl)) return(FALSE)
  ok <- tryCatch(parallel::clusterCall(cl, function() TRUE),
                 error = function(e) NULL)
  !is.null(ok) && length(ok) == length(cl) &&
    all(vapply(ok, isTRUE, logical(1)))
}

.par_stop_cluster <- function() {
  cl <- .par_state$cl
  if (!is.null(cl)) tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
  .par_state$cl <- NULL
  .par_state$n  <- NULL
  invisible(NULL)
}

# A ready pool of exactly `n` workers with the package loaded, or NULL if one
# could not be built. Reuses the cached pool when it is the right size and
# still answering; rebuilds otherwise, since a pool of the wrong size would
# silently change how the work is divided.
.par_cluster <- function(n) {
  if (identical(.par_state$n, n) && .par_cluster_alive(.par_state$cl))
    return(.par_state$cl)
  .par_stop_cluster()

  cl <- tryCatch(parallel::makePSOCKcluster(n), error = function(e) NULL)
  if (is.null(cl)) return(NULL)

  dev <- .par_dev_path()
  ok <- tryCatch({
    if (is.null(dev)) {
      parallel::clusterCall(cl, function(lp) {
        .libPaths(lp)
        requireNamespace("mixtureEM", quietly = TRUE)
      }, .libPaths())
    } else {
      parallel::clusterCall(cl, function(lp, p) {
        .libPaths(lp)
        requireNamespace("pkgload", quietly = TRUE) &&
          isTRUE(tryCatch({ pkgload::load_all(p, quiet = TRUE); TRUE },
                          error = function(e) FALSE))
      }, .libPaths(), dev)
    }
  }, error = function(e) list(FALSE))

  if (!all(vapply(ok, isTRUE, logical(1)))) {
    tryCatch(parallel::stopCluster(cl), error = function(e) NULL)
    return(NULL)
  }

  .par_state$cl <- cl
  .par_state$n  <- n
  cl
}

.onUnload <- function(libpath) .par_stop_cluster()

# lapply() over `X`, on `n_cores` worker processes when that is worth doing.
#
# PSOCK rather than forking: the package is used heavily on Windows, which has
# no fork(), and a PSOCK cluster behaves the same everywhere. The cost is that
# each worker is a fresh session that must load the package and be sent its
# share of the data. The session is paid for once (see .par_cluster() above);
# the data is sent on every call, which is why this stays opt-in and remains
# pointless for a job whose pieces are shorter than the transfer.
.par_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores)) n_cores <- 1L
  if (n_cores <= 1L || length(X) < 2L) return(lapply(X, FUN, ...))

  # Sized by what the caller asked for, never by how many pieces this
  # particular call has. The pool is cached across calls, and a pool rebuilt
  # every time the job length changed would pay the spawn cost it exists to
  # avoid -- a five-restart fit followed by a twenty-restart one would rebuild
  # twice. Workers with nothing to do cost memory and no time.
  n_cores <- min(n_cores, .par_max_cores())
  if (n_cores <= 1L) return(lapply(X, FUN, ...))

  cl <- .par_cluster(n_cores)
  if (is.null(cl)) {
    warning("could not start ", n_cores, " workers; running sequentially",
            call. = FALSE)
    return(lapply(X, FUN, ...))
  }

  # parLapplyLB() rather than parLapply(): the latter cuts `X` into one fixed
  # chunk per worker before any of it runs, which is only efficient when the
  # pieces cost the same. The pieces here are model fits, and they do not: two
  # random starts of the same model routinely differ several-fold in how many
  # iterations they need, and a bootstrap replicate is no more predictable.
  # Under static chunking the call cannot finish before its unluckiest chunk,
  # so the workers that drew the cheap pieces wait. parLapplyLB() hands out one
  # piece at a time to whichever worker is free, so the tail is one piece long
  # instead of one chunk long.
  #
  # This is worth having wherever the pool is smaller than the job -- a BLRT's
  # two thousand replicates, a hundred-restart search. It cannot help a call
  # with fewer pieces than workers, which is a different problem with a
  # different fix.
  #
  # Nothing about the answer changes. Which worker runs a restart was never an
  # input to it -- workers draw no random numbers, see the note at the top of
  # this file -- and clusterApplyLB() returns results in the order of `X`, not
  # the order they finished. The bit-for-bit equality across `n_cores` that the
  # validation targets rest on is untouched.
  #
  # The pool outlives this call, so a failure here must not leave it holding a
  # half-finished state. stopCluster() on error is the safe reading: the next
  # call rebuilds, which costs a spawn and cannot cost a wrong answer.
  res <- tryCatch(parallel::parLapplyLB(cl, X, FUN, ...),
                  error = function(e) { .par_stop_cluster(); stop(e) })
  res
}
