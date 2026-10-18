# ==============================================================================
# The one-score class: the package's single definition of a degenerate fit
# ==============================================================================
#
# A mixture of normals whose classes carry variances of their own has an
# unbounded likelihood: let one class own a set of cases that share one exact
# value and drive that class's variance towards zero, and the density at those
# cases grows without limit (Kiefer & Wolfowitz, 1956; Day, 1969). Real data
# hand the likelihood such sets all the time -- a scale's floor or ceiling, a
# rounded mean score, a column of zeros -- and an EM start that walks into one
# outscores every meaningful solution.
#
# This file states that condition once, in data terms, for every family that
# can meet it: a class whose own variance governs a set of observations of which
# more than `.one_score_share` of the class's posterior mass sits on one exact
# value. "Its own" is literal. A variance held equal across the classes is
# never tested (the likelihood is bounded in that direction, and a class sitting
# on a pile of tied scores is a legitimate solution there); a variance shared
# across occasions, groups or blocks is tested on everything it governs, pooled,
# because one occasion on a single value cannot send a variance that the other
# occasions also inform to zero.
#
# The definition does not depend on the prior. Under the default variance prior
# the collapsing variance stops at a finite floor instead of at zero, but the
# class still sits on one value; with the prior off it runs to the floor. So
# one rule serves both, which a threshold on the variance itself cannot: the
# collapsed variance under the prior (about 1e-3 of the marginal) overlaps the
# genuine small classes of a priors-off fit.
#
# It is what every family's search sets aside (.em_state_degenerate(),
# fit_lta()'s ranking) and what the collapse warnings report
# (.check_gaussian_degeneracy(), .gmm_boundary()). The 90% cut sits in a wide
# empty gap: on simulated profile and growth data with floor and ceiling
# effects, and on two published latent profile data sets, every admissible
# start measured put at most 68% of a class on one value and every collapsed
# one at least 99.8%.

.one_score_share <- 0.9

# The class-specific variance parameters of an emission, each as the class and
# the (column, block) pairs of the data it governs. `offset` and `block` place
# a sub-model inside a blocks or nested model. Families with no class-specific
# variance return an empty list.
.dispersion_groups <- function(mm, offset = 0L, block = 1L) {
  if (is.null(mm)) return(list())
  K <- mm$n_components
  group <- function(k, cols, blocks = rep(block, length(cols)))
    list(class = k, cols = as.integer(cols), blocks = as.integer(blocks))

  if (inherits(mm, "blocks")) {
    J  <- mm$n_items
    gs <- unlist(lapply(seq_len(mm$n_blocks), function(b)
      .dispersion_groups(mm$models[[b]], offset + (b - 1L) * J, b)),
      recursive = FALSE)
    if (!length(gs)) return(list())
    # An item whose variance the blocks share is one parameter per class: pool
    # its cells across the blocks.
    shared <- if ("covariances" %in% mm$invariant_params) seq_len(J)
              else mm$invariant_items
    if (!length(shared)) return(gs)
    item <- vapply(gs, function(g) (g$cols[1] - offset - 1L) %% J + 1L, 1L)
    cls  <- vapply(gs, `[[`, 1L, "class")
    key  <- ifelse(item %in% shared, paste(cls, item), paste0("u", seq_along(gs)))
    return(lapply(split(gs, factor(key, unique(key))), function(s)
      group(s[[1]]$class, unlist(lapply(s, `[[`, "cols")),
            unlist(lapply(s, `[[`, "blocks")))))
  }

  if (inherits(mm, "nested")) {
    off <- offset
    out <- list()
    for (nm in names(mm$models)) {
      out <- c(out, .dispersion_groups(mm$models[[nm]], off, block))
      off <- off + as.integer(mm$columns_per_model[[nm]])
    }
    return(out)
  }

  if (class(mm)[1] %in% c("gaussian_diag", "gaussian_diag_nan")) {
    if (isTRUE(mm$variances_equal)) return(list())
    J <- ncol(mm$parameters$means %||% matrix(0, 0, 0))
    return(unlist(lapply(seq_len(K), function(k) lapply(seq_len(J), function(j)
      group(k, offset + j))), recursive = FALSE))
  }

  if (inherits(mm, "structured_normal")) {
    psi_k   <- !isTRUE(mm$psi_equal) && length(mm$r_cols) > 0L
    theta_k <- !isTRUE(mm$theta_equal)
    if (!psi_k && !theta_k) return(list())
    Tn <- nrow(mm$design)
    # Only a class-specific residual variance constant over occasions makes the
    # class's dispersion one parameter for all of them.
    pooled <- !psi_k && isTRUE(mm$theta_shared_occasions)
    return(unlist(lapply(seq_len(K), function(k)
      if (pooled) list(group(k, offset + seq_len(Tn)))
      else lapply(seq_len(Tn), function(t) group(k, offset + t))),
      recursive = FALSE))
  }

  if (inherits(mm, "lcga") && isTRUE(mm$fam$has_dispersion)) {
    Tn <- nrow(mm$design)
    return(lapply(seq_len(K), function(k) group(k, offset + seq_len(Tn))))
  }

  list()
}

# The groups of `mm` whose class sits on one value of the data, as a data frame
# with the class, the column holding that value, the value and the share; NULL
# when there are none. `post` is the n x K posterior matrix, or a list of them
# indexed by block (a transition model's status posteriors, one per occasion).
.one_score_cells <- function(mm, X, post, weights = NULL,
                             share = .one_score_share) {
  gs <- .dispersion_groups(mm)
  if (!length(gs) || is.null(post) || is.null(X)) return(NULL)
  n <- nrow(X)
  w <- if (is.null(weights) || length(weights) != n) rep(1, n) else weights
  post_of <- function(b) if (is.list(post)) post[[b]] else post

  rows <- lapply(gs, function(g) {
    vals <- numeric(0); mass <- numeric(0); where <- integer(0)
    for (i in seq_along(g$cols)) {
      j <- g$cols[i]
      if (j > ncol(X)) next
      P <- post_of(g$blocks[i])
      if (is.null(P) || g$class > ncol(P)) next
      ok <- !is.na(X[, j])
      vals  <- c(vals, X[ok, j])
      mass  <- c(mass, w[ok] * P[ok, g$class])
      where <- c(where, rep(j, sum(ok)))
    }
    total <- sum(mass)
    if (!length(vals) || !is.finite(total) || total <= 0) return(NULL)
    u   <- unique(vals)
    idx <- match(vals, u)
    by_value <- rowsum(mass, idx)
    top <- which.max(by_value)
    s <- by_value[top] / total
    if (!is.finite(s) || s <= share) return(NULL)
    t_val <- as.integer(rownames(by_value)[top])
    data.frame(class = g$class, col = where[match(t_val, idx)], value = u[t_val],
               share = s)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out[order(-out$share), , drop = FALSE]
}
