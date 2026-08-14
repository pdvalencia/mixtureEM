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
  cat("  These estimates are not interpretable, and this fit's BIC cannot be\n")
  cat("  compared with a clean one's. See ?fit_mixture for what to do.\n\n")
  invisible(NULL)
}

# The number of classes the fit estimated, used to state the prior remedy in the
# user's own K rather than as a bare constant. `n_components` on a
# mixture_model, `n_statuses` on an lta_model.
.degeneracy_n_classes <- function(fit) {
  k <- fit$n_components %||% fit$n_statuses
  if (is.numeric(k) && length(k) == 1L && is.finite(k) && k >= 1) as.integer(k)
  else NA_integer_
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

  # The prior remedy is stated as one artificial observation per class, and
  # printed as the number that means for *this* model. A bare constant does not
  # transfer: the constant is spread over the classes, so the same value is a
  # different amount of prior at every K. Calibrated against a reference
  # implementation on five-point scales, one observation per class was the
  # weakest setting that both lifted the flagged variance into the range of the
  # model's genuinely small variances and moved its class mean off the scale
  # ceiling. See the `bayes_constants` section of ?fit_mixture.
  K <- .degeneracy_n_classes(fit)
  prior_hint <- if (is.na(K)) "bayes_constants = list(variances = <n_classes>)"
                else sprintf("bayes_constants = list(variances = %d)", K)

  # R truncates a condition message at getOption("warning.length"), which
  # defaults to 1000 bytes and which RStudio does not raise -- so a warning
  # written to full length loses its tail in the console most applied users are
  # in, and the tail is where the remedies are. This message is therefore kept
  # under that limit at the worst case .gaussian_boundary_lines() can produce,
  # and the reasoning it used to carry lives in ?fit_mixture instead. The
  # asymmetry with .print_degenerate_note() is deliberate: cat() has no limit, so
  # the printed note and the help file are the long form and this is the short
  # one. The flagged cells come first because they are the only content the user
  # cannot get anywhere else, and so must never be what gets cut.
  warning(sprintf(
    paste0("A class variance has collapsed towards zero: %s. ",
           "These estimates are not interpretable, and this fit's BIC cannot ",
           "be compared with a clean fit's. Three ways out, to choose between ",
           "on substantive grounds: (1) variances_equal = TRUE; (2) fewer ",
           "classes; or (3) a stronger prior, %s. ",
           "See ?fit_mixture for why, and what to check afterwards."),
    paste(.gaussian_boundary_lines(flagged), collapse = "; "), prior_hint),
    call. = FALSE)
  fit
}
