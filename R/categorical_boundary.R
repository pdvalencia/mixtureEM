# ==============================================================================
# Boundary detection for categorical item probabilities under a continuous
# random intercept
# ==============================================================================
#
# The Gaussian collapse in R/gaussian_boundary.R has a categorical counterpart,
# but it lives in a different place than the obvious one. A plain categorical
# M-step (m_step.bernoulli / m_step.multinoulli, R/categorical.R) estimates each
# class-by-item probability as a weighted count ratio, bounded in [0, 1] by
# construction; with no prior mass it can sit arbitrarily close to 0 or 1 as a
# stable, replicated, externally-matched optimum whenever a response never
# co-occurs with a class in the weighted sample -- there is no logit-scale
# ridge for it to run away on. Measured directly (RECORDS.md, "R11" section 8):
# a known-good, externally-validated fit reaches a cell probability of 1e-7 on
# 2587 cases, so a threshold on mm$parameters$pis would flag fits already known
# to be correct. This file therefore checks only the continuous random
# intercept's item parameters (fit$ri$A for binary indicators, fit$ri$theta for
# ordinal), which are logit-scale intercepts run through a numerically
# integrated random factor with nothing to bound them directly -- exactly the
# mechanism that produced the -14441.9787 degenerate LTA-FAQ optimum this guard
# exists for (RECORDS.md, "R11" sections 6-7).

# One row per class/item (binary) or class/threshold-column (ordinal). For
# ordinal indicators every column of ri$theta is checked on the same |value| >
# threshold rule: the first column per item is a logit-scale threshold and the
# remaining ones are log-scale gaps to the next threshold (the increment
# parameterisation, R/ordinal.R), so a large value in either means the same
# thing -- a category probability collapsing towards 0.
.ri_logit_cells <- function(ri, item_names = NULL) {
  if (is.null(ri)) return(NULL)
  if (!is.null(ri$theta)) {
    M <- ri$theta
    cats <- ri$cats
    col_item <- if (!is.null(cats))
      rep(seq_along(cats), times = cats - 1L) else rep(1L, ncol(M))
  } else if (!is.null(ri$A)) {
    M <- ri$A
    col_item <- seq_len(ncol(M))
  } else return(NULL)

  nms <- item_names %||% paste0("item ", seq_len(max(col_item)))
  data.frame(class = as.vector(row(M)), col = as.vector(col(M)),
             item  = nms[col_item[as.vector(col(M))]],
             logit = as.vector(M))
}

# Flag continuous-RI item-parameter cells past the boundary threshold.
#
# threshold is stated on |logit|: a known-clean continuous-RI fit tops out at
# 4.301, a known-degenerate one starts at 17.976 and grows without bound under
# continued EM (RECORDS.md, "R11" section 8) -- 8 sits with an order of
# magnitude of slack on both sides and is also where the problem first became
# visible by inspection.
#
# Returns NULL when nothing is flagged, so callers can test with length() --
# in particular, NULL for every fit with no random intercept. `X` is accepted
# and ignored, only to keep the same call shape as .gaussian_boundary(fit$mm,
# X, ...) at the one call site both share (.check_gaussian_degeneracy()).
.categorical_boundary <- function(fit, X = NULL, threshold = 8) {
  item_names <- fit$longitudinal$item_names
  cells <- .ri_logit_cells(fit$ri, item_names)
  if (is.null(cells) || !nrow(cells)) return(NULL)

  flagged <- cells[is.finite(cells$logit) & abs(cells$logit) > threshold, ,
                   drop = FALSE]
  if (!nrow(flagged)) return(NULL)

  flagged$prob <- plogis(flagged$logit)
  flagged$kind <- "probability"
  flagged[order(-abs(flagged$logit)),
          c("class", "col", "item", "logit", "prob", "kind")]
}

# Format the flagged cells for a warning or for print(). Mirrors
# .gaussian_boundary_lines().
.categorical_boundary_lines <- function(flagged, max_show = 5L) {
  n <- nrow(flagged)
  show <- flagged[seq_len(min(n, max_show)), , drop = FALSE]
  lines <- sprintf("%s in class %d (logit %+.2f, p = %.3g)",
                   show$item, show$class, show$logit, show$prob)
  if (n > max_show)
    lines <- c(lines, sprintf("and %d more", n - max_show))
  lines
}
