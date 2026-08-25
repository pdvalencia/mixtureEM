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

# lapply() over `X`, on `n_cores` worker processes when that is worth doing.
#
# PSOCK rather than forking: the package is used heavily on Windows, which has
# no fork(), and a PSOCK cluster behaves the same everywhere. The cost is that
# each worker is a fresh session that must load the package and be sent its
# share of the data, which is why this is opt-in and pointless for short jobs.
.par_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores)) n_cores <- 1L
  if (n_cores <= 1L || length(X) < 2L) return(lapply(X, FUN, ...))

  n_cores <- min(n_cores, length(X), .par_max_cores())
  if (n_cores <= 1L) return(lapply(X, FUN, ...))

  cl <- tryCatch(parallel::makePSOCKcluster(n_cores), error = function(e) NULL)
  if (is.null(cl)) {
    warning("could not start ", n_cores, " workers; running sequentially",
            call. = FALSE)
    return(lapply(X, FUN, ...))
  }
  on.exit(parallel::stopCluster(cl), add = TRUE)

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
    warning("workers could not load the package; running sequentially",
            call. = FALSE)
    return(lapply(X, FUN, ...))
  }

  parallel::parLapply(cl, X, FUN, ...)
}
