# ==============================================================================
# Collapsed-variance detection for continuous indicators
# ==============================================================================
#
# A mixture of normals with free variances has an unbounded likelihood: drive one
# class's variance on one item towards zero with its mean sitting on a handful of
# near-identical cases, and the density at those cases grows without limit. The
# prior in m_step.gaussian_diag() and in refine_lbfgs() bounds the objective, so
# the estimator no longer runs off to infinity — but at the default prior
# strength the collapsed path remains a *finite* spurious optimum, and a spurious
# optimum that wins on the objective is reported like any other solution.
#
# That is the reason detection is not optional. The choice made here is to
# default to a weak prior, whose failure mode is loud and recognisable, rather
# than a strong one, which would suppress this particular failure at the cost of
# quietly shifting perfectly healthy small-sample fits. A default that fails
# loudly is worth more than one that fails silently — but only if something is
# listening, which is what this file is.
#
# The signature is unmistakable once looked for: a variance orders of magnitude
# below the item's own marginal, usually with the class mean pinned exactly at
# the smallest or largest observed value of that item.

# Every class-by-item variance in a measurement model, mapped back to the data
# column it belongs to.
#
# Walks the three shapes a Gaussian emission can arrive in: flat, a `blocks`
# model (repeated measures or multiple groups, one sub-model per block), and a
# `nested` model (one sub-model per column range, of which only some are
# continuous). Anything that is not a `gaussian_diag` family contributes
# nothing, which is also how a mixed measurement model gets handled: its binary
# and count blocks are simply skipped.
.gaussian_variance_cells <- function(mm, offset = 0L) {
  if (is.null(mm)) return(NULL)

  if (inherits(mm, "blocks")) {
    J <- mm$n_items
    return(do.call(rbind, lapply(seq_len(mm$n_blocks), function(b)
      .gaussian_variance_cells(mm$models[[b]], offset + (b - 1L) * J))))
  }

  if (inherits(mm, "nested")) {
    off <- offset
    out <- vector("list", length(mm$models))
    for (i in seq_along(mm$models)) {
      out[[i]] <- .gaussian_variance_cells(mm$models[[i]], off)
      off <- off + as.integer(mm$columns_per_model[[i]])
    }
    return(do.call(rbind, out))
  }

  if (!class(mm)[1] %in% c("gaussian_diag", "gaussian_diag_nan")) return(NULL)
  V <- mm$parameters$covariances
  M <- mm$parameters$means
  if (is.null(V) || !length(V)) return(NULL)

  data.frame(class    = as.vector(row(V)),
             col      = offset + as.vector(col(V)),
             variance = as.vector(V),
             mean     = as.vector(M))
}

# Flag class/item cells whose variance has collapsed.
#
# The threshold is a *ratio* to the item's observed marginal variance, not an
# absolute value, so that it means the same thing on a five-point scale and on an
# income variable. At 1% it has an order of magnitude of slack on both sides of
# everything measured while this was being diagnosed: genuinely small but healthy
# class variances ran 2-7% of the marginal, and collapsed ones 0.03-0.08%.
#
# Returns NULL when nothing is flagged, so callers can test with length().
.gaussian_boundary <- function(mm, X, threshold = 0.01, weights = NULL) {
  cells <- .gaussian_variance_cells(mm)
  if (is.null(cells) || !nrow(cells)) return(NULL)

  cells <- cells[cells$col >= 1L & cells$col <= ncol(X), , drop = FALSE]
  if (!nrow(cells)) return(NULL)

  # The same weighted marginal the prior is centred on. Under a survey design an
  # unweighted marginal would put the threshold on a different scale from the
  # quantity being compared to it.
  if (!is.null(weights) && length(weights) != nrow(X)) weights <- NULL
  s2 <- .marginal_var(X, weights)
  cells$marginal <- s2[cells$col]
  cells$ratio    <- cells$variance / cells$marginal

  flagged <- cells[is.finite(cells$ratio) & cells$ratio < threshold, ,
                   drop = FALSE]
  if (!nrow(flagged)) return(NULL)

  # Corroborating signal: the class mean sitting on the smallest or largest
  # observed value of the item. A collapsed class is usually a clump of cases at
  # a scale endpoint, so this is what distinguishes "the model found a floor
  # effect and called it a class" from an ordinary tight class. Reported when it
  # holds rather than required, since a collapse in the middle of a scale is
  # equally spurious.
  rng <- apply(X, 2, range, na.rm = TRUE)
  tol <- pmax((rng[2, ] - rng[1, ]) * 1e-4, 1e-8)
  flagged$pinned <-
    abs(flagged$mean - rng[1, flagged$col]) <= tol[flagged$col] |
    abs(flagged$mean - rng[2, flagged$col]) <= tol[flagged$col]

  nms <- colnames(X) %||% paste0("column ", seq_len(ncol(X)))
  flagged$item <- nms[flagged$col]

  flagged[order(flagged$ratio), c("class", "col", "item", "variance",
                                  "marginal", "ratio", "mean", "pinned")]
}

# Format the flagged cells for a warning or for print().
.gaussian_boundary_lines <- function(flagged, max_show = 5L) {
  n <- nrow(flagged)
  show <- flagged[seq_len(min(n, max_show)), , drop = FALSE]
  lines <- sprintf("%s in class %d (variance %.3g vs %.3g for the item overall%s)",
                   show$item, show$class, show$variance, show$marginal,
                   ifelse(show$pinned, "; class mean at a data boundary", ""))
  if (n > max_show)
    lines <- c(lines, sprintf("and %d more", n - max_show))
  lines
}

# Repeat the flag in print()/summary(). No-op for a clean fit.
.print_degenerate_note <- function(x) {
  flagged <- x$degenerate
  if (is.null(flagged) || !nrow(flagged)) return(invisible(NULL))
  cat("\nWARNING - collapsed class variance:\n")
  for (line in .gaussian_boundary_lines(flagged)) cat("  ", line, "\n", sep = "")
  cat("  These estimates are not interpretable. Refit with\n")
  cat("  bayes_constants = list(variances = 5), or with fewer classes.\n\n")
  invisible(NULL)
}

# Run the check on a fitted model, store the result, and warn.
#
# Called from both fitting paths -- fit_mixture_internal() and fit_lta(), which
# has its own EM driver -- so a continuous indicator is checked wherever it is
# estimated. Silent for every model with no continuous indicators.
.check_gaussian_degeneracy <- function(fit, X, quiet = FALSE) {
  # `sample_weights` on a mixture_model, `weights_vec` on an lta_model.
  w <- fit$sample_weights %||% fit$weights_vec
  flagged <- tryCatch(.gaussian_boundary(fit$mm, X, weights = w),
                      error = function(e) NULL)
  fit$degenerate <- flagged
  if (is.null(flagged) || isTRUE(quiet)) return(fit)

  warning(sprintf(
    paste0("A class variance has collapsed towards zero: %s. ",
           "The likelihood of a mixture of normals is unbounded in this ",
           "direction, so this solution can score better than any meaningful ",
           "one while describing a handful of near-identical cases rather ",
           "than a subgroup. Do not interpret it as it stands. Remedies, in ",
           "the order worth trying: refit with ",
           "bayes_constants = list(variances = 5), which strengthens the ",
           "prior holding variances away from zero; fit fewer classes; or ",
           "inspect the distribution of the named item for a floor, ceiling ",
           "or spike that a class has latched onto."),
    paste(.gaussian_boundary_lines(flagged), collapse = "; ")),
    call. = FALSE)
  fit
}
