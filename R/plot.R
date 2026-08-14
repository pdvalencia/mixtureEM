# ==============================================================================
# Profile plot for fitted mixture models (base R, zero dependencies)
# ==============================================================================
# A single-panel "profile plot" places every indicator on a common [0, 1] axis
# so that classes can be compared at a glance:
#
#   * Binary indicators      -> probability of endorsement (already 0-1).
#   * Continuous indicators  -> class means, min-max scaled against the
#                               observed data range.
#   * Polytomous indicators  -> expected (mean) category, min-max scaled to
#                               [0, 1] over the 1..M category range.
#
# Scaled items are flagged with a trailing "*" and a footnote, since their
# position reflects relative standing rather than a probability.

# Okabe-Ito qualitative palette: the standard choice for colour-vision
# deficiency. Recycled if a model has more than eight classes.
.okabe_ito <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermilion
  "#CC79A7", # reddish purple
  "#000000"  # black
)

# Legend labels for class trajectories/profiles: the user's label (capped so a
# verbose one cannot push the legend off the device) plus the class share.
# Shared by every plot method that draws one line per class.
.class_plot_labels <- function(base_labels, weights, max_chars = 24L) {
  base <- .shorten_labels(as.character(base_labels), width = max_chars)
  sprintf("%s (%.1f%%)", as.vector(base), weights * 100)
}

# Does a (possibly nested) measurement model contain any block plotted on the
# observed data range? Continuous means and Poisson rates are both rescaled
# against the raw indicators, so both need the raw data to be available.
.has_continuous <- function(mm) {
  if (inherits(mm, c("nested", "blocks")))
    return(any(vapply(mm$models, .has_continuous, logical(1))))
  !is.null(mm$parameters$means) || !is.null(mm$parameters$rates)
}

# Build a [classes x items] matrix of values on a common [0, 1] scale.
# Column names carry a trailing "*" for any item that was rescaled. For nested
# (mixed-type) models the blocks are concatenated in submodel order.
.prepare_profile_data <- function(mm, indicators = NULL) {

  if (inherits(mm, c("nested", "blocks"))) {
    blocks <- lapply(mm$models, .prepare_profile_data, indicators = indicators)
    return(do.call(cbind, blocks))
  }

  is_cont  <- !is.null(mm$parameters$means)
  is_poly  <- !is.null(mm$max_val)
  is_bin   <- !is.null(mm$parameters$pis) && is.null(mm$max_val)
  is_count <- !is.null(mm$parameters$rates)

  plot_data  <- list()

  # --- Continuous: min-max scale class means to [0, 1] ----------------------
  if (is_cont) {
    means <- mm$parameters$means
    cols  <- colnames(means) %||% paste0("Cont_", seq_len(ncol(means)))

    for (j in seq_len(ncol(means))) {
      # Prefer the observed data range; fall back to the range of the class
      # means only when the raw data is unavailable.
      have_obs <- !is.null(indicators) && cols[j] %in% colnames(indicators)
      if (have_obs) {
        col_min <- min(indicators[, cols[j]], na.rm = TRUE)
        col_max <- max(indicators[, cols[j]], na.rm = TRUE)
      } else {
        col_min <- min(means[, j]); col_max <- max(means[, j])
      }
      denom <- if (col_max == col_min) 1 else (col_max - col_min)
      plot_data[[paste0(cols[j], "*")]] <- (means[, j] - col_min) / denom
    }
  }

  # --- Binary: probabilities are already on [0, 1] --------------------------
  if (is_bin) {
    pis  <- mm$parameters$pis
    cols <- colnames(pis) %||% paste0("Bin_", seq_len(ncol(pis)))
    for (j in seq_len(ncol(pis)))
      plot_data[[cols[j]]] <- pis[, j]
  }

  # --- Count: rates share the continuous treatment (min-max to [0, 1]) ------
  # Rates live on the data's own scale, like continuous means, so the same
  # rescaling puts them on the common profile axis.
  if (is_count) {
    rates <- mm$parameters$rates
    cols  <- colnames(rates) %||% paste0("Count_", seq_len(ncol(rates)))

    for (j in seq_len(ncol(rates))) {
      have_obs <- !is.null(indicators) && cols[j] %in% colnames(indicators)
      if (have_obs) {
        col_min <- min(indicators[, cols[j]], na.rm = TRUE)
        col_max <- max(indicators[, cols[j]], na.rm = TRUE)
      } else {
        col_min <- min(rates[, j]); col_max <- max(rates[, j])
      }
      denom <- if (col_max == col_min) 1 else (col_max - col_min)
      plot_data[[paste0(cols[j], "*")]] <- (rates[, j] - col_min) / denom
    }
  }

  # --- Polytomous (ordinal): expected category, scaled to [0, 1] ------------
  if (is_poly) {
    pis     <- mm$parameters$pis
    M       <- mm$max_val
    n_items <- ncol(pis) / M
    base    <- mm$item_names %||% paste0("Ord_", seq_len(n_items))
    if (length(base) != n_items) base <- paste0("Ord_", seq_len(n_items))

    categories <- seq_len(M)
    for (j in seq_len(n_items)) {
      item_probs <- pis[, ((j - 1) * M + 1):(j * M), drop = FALSE]
      expected   <- rowSums(sweep(item_probs, 2, categories, "*"))
      # Expected value ranges over [1, M]; map to [0, 1].
      plot_data[[paste0(base[j], "*")]] <- (expected - 1) / (M - 1)
    }
  }

  if (length(plot_data) == 0L) return(NULL)
  do.call(cbind, plot_data)
}

#' Profile Plot for a Fitted Mixture Model
#'
#' @description
#' Draws a profile plot of the measurement model: one line per latent class,
#' with every indicator placed on a common \[0, 1\] axis. Binary indicators are
#' shown as endorsement probabilities; continuous indicators are min-max scaled
#' against their observed range; polytomous indicators are summarised by their
#' expected category and scaled to \[0, 1\]. Rescaled items are marked with "*".
#'
#' Uses only base graphics and the colour-blind-friendly Okabe-Ito palette.
#'
#' @param x A fitted \code{mixture_model} object.
#' @param main Plot title.
#' @param class_labels Optional character vector of labels for the classes.
#'   Defaults to "Class 1", "Class 2", ...
#' @param colors Optional vector of colours (one per class). Defaults to the
#'   Okabe-Ito palette, recycled if necessary.
#' @param ... Currently unused; present for S3 compatibility.
#'
#' @return The fitted model, invisibly.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
#' plot(fit)
#'
#' @importFrom graphics par matplot axis text legend mtext
#' @export
plot.mixture_model <- function(x, main = "Latent Class / Profile Plot",
                               class_labels = NULL, colors = NULL, ...) {

  if (is.null(x$mm)) stop("No measurement model to plot.")

  if (is.null(x$data) && .has_continuous(x$mm))
    message("Raw indicators are not stored on this fit, so continuous items ",
            "are scaled against their estimated class means rather than the ",
            "observed data range. Refit with the current version to enable ",
            "observed-range scaling.")

  plot_mat <- .prepare_profile_data(x$mm, indicators = x$data)
  if (is.null(plot_mat) || ncol(plot_mat) == 0L)
    stop("Could not extract plottable parameters from the measurement model.")

  n_classes <- nrow(plot_mat)
  n_items   <- ncol(plot_mat)

  my_colors <- if (is.null(colors)) rep(.okabe_ito, length.out = n_classes)
               else rep(colors, length.out = n_classes)
  my_shapes <- rep(15:20, length.out = n_classes)

  base_labels <- if (is.null(class_labels)) paste("Class", seq_len(n_classes))
                 else class_labels

  labels <- .class_plot_labels(base_labels, x$weights)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  # Margins sized from the labels that actually go in them: the right margin
  # holds the legend (whose labels are already capped by .class_plot_labels),
  # the bottom margin the 45-degree item names. A fixed margin either clips
  # long labels or wastes space on short ones.
  mar_right  <- min(3 + 0.45 * max(nchar(labels)), 18)
  mar_bottom <- max(6, min(3.5 + 0.4 * max(nchar(colnames(plot_mat))), 12))
  # The first item's rotated label leans left of the axis, so the left margin
  # grows with it too (the y-axis title needs the 4-line floor regardless).
  mar_left   <- max(4, min(2 + 0.25 * nchar(colnames(plot_mat)[1]), 8))
  par(mar = c(mar_bottom, mar_left, 4, mar_right), xpd = TRUE)

  matplot(
    x    = seq_len(n_items),
    y    = t(plot_mat),
    type = "b",
    pch  = my_shapes,
    lty  = 1,
    lwd  = 2,
    col  = my_colors,
    ylim = c(0, 1),
    xlim = c(0.75, n_items + 0.25),
    xaxt = "n",
    xlab = "",
    ylab = "Scaled value / probability",
    main = main,
    bty  = "l",
    las  = 1
  )

  axis(1, at = seq_len(n_items), labels = FALSE)
  text(
    x      = seq_len(n_items),
    y      = par("usr")[3] - 0.05,
    labels = colnames(plot_mat),
    srt    = 45,
    adj    = 1,
    cex    = 0.9
  )

  # Offset in a fixed fraction of the x-range, not a fixed data-unit amount:
  # 0.1 data units is a different physical distance for a 4-item plot than for
  # a 30-item one, which used to push the legend off the device.
  x_pos <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  y_pos <- par("usr")[4]

  legend(
    x      = x_pos,
    y      = y_pos,
    legend = labels,
    col    = my_colors,
    pch    = my_shapes,
    lty    = 1,
    lwd    = 2,
    bty    = "n",
    cex    = 0.9,
    title  = "Class"
  )

  if (any(grepl("\\*", colnames(plot_mat))))
    mtext(
      "* Continuous items min-max scaled; ordinal items shown as scaled expected category.",
      side = 1, line = mar_bottom - 1.2, adj = 0, cex = 0.8, font = 3,
      col = "grey30"
    )

  invisible(x)
}

# ==============================================================================
# Elbow plot for a compare_mixtures() / compare_longitudinal() sweep
# ==============================================================================

#' Elbow Plot for a Model-Selection Sweep
#'
#' @description
#' Plots the information criteria from [`compare_mixtures()`] or
#' [`compare_longitudinal()`] against the number of classes, which is the
#' picture the fit table is usually read as. Masyn (2013) describes the reading:
#' plot the criteria against K and look for the point of diminishing returns
#' rather than taking the raw minimum, because the BIC often keeps falling
#' slowly as classes are added without those classes being substantively
#' distinct. The minimum is marked so it can be seen, not so it can be obeyed.
#'
#' The plot informs the decision; it does not make it. Ram and Grimm (2009,
#' p. 571) put it that "model selection is an art" — the criteria are one input
#' alongside class size, interpretability and the substantive question.
#'
#' @param x A `mixture_comparison` object, as returned by
#'   [`compare_mixtures()`] or [`compare_longitudinal()`].
#' @param indices Character vector, any subset of `"BIC"`, `"AIC"` and
#'   `"SABIC"`. Default `"BIC"` alone. All three are \eqn{-2\ell} plus a penalty
#'   and so share one axis; the log-likelihood and the entropy are deliberately
#'   not allowed on it, being on a different scale entirely.
#' @param entropy Logical. When `TRUE`, relative entropy is drawn in a second
#'   panel below, on a fixed 0-1 axis. It gets its own panel rather than a
#'   right-hand axis, since a twin axis would invite reading the two against
#'   each other, which is exactly the comparison it must not support: entropy
#'   measures how well separated the classes are, not how well the model fits,
#'   and it is not a model-selection criterion.
#' @param main Optional title for the top panel.
#' @param ... Currently unused. Present for S3 method compatibility.
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @references
#' Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. In
#' T. D. Little (Ed.), \emph{The Oxford Handbook of Quantitative Methods}
#' (Vol. 2, pp. 551-611). Oxford University Press.
#'
#' Ram, N., & Grimm, K. J. (2009). Growth mixture modeling: A method for
#' identifying differences in longitudinal change among unobserved groups.
#' \emph{International Journal of Behavioral Development}, \emph{33}(6),
#' 565-576. \doi{10.1177/0165025409343765}
#'
#' @seealso [`compare_mixtures()`], [`compare_longitudinal()`]
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
#' result <- compare_mixtures(X, k_range = 1:4, measurement = "binary",
#'                            n_init = 5)
#' plot(result)
#' plot(result, indices = c("BIC", "AIC"), entropy = TRUE)
#'
#' @export
plot.mixture_comparison <- function(x, indices = "BIC", entropy = FALSE,
                                    main = NULL, ...) {
  tab <- x$fit_table
  if (is.null(tab) || !nrow(tab))
    stop("Nothing to plot: the comparison has no fit table.", call. = FALSE)

  indices <- match.arg(indices, c("BIC", "AIC", "SABIC"), several.ok = TRUE)
  missing_cols <- setdiff(indices, names(tab))
  if (length(missing_cols))
    stop(sprintf("The fit table has no %s column.",
                 paste(missing_cols, collapse = " or ")), call. = FALSE)

  K     <- tab$Classes
  # `%in% TRUE` rather than the column itself, so an NA cell reads as "not
  # flagged" instead of propagating into every pch it touches.
  unrep <- if (is.null(tab$Unreplicated)) rep(FALSE, length(K))
           else tab$Unreplicated %in% TRUE

  # Colour is never the only channel that distinguishes the lines: line type and
  # plotting symbol vary with it, so the plot survives a greyscale print and a
  # colour-vision-deficient reader.
  cols     <- rep(.okabe_ito, length.out = length(indices))
  ltys     <- rep(c(1, 2, 4),   length.out = length(indices))
  pch_fill <- rep(c(16, 17, 15), length.out = length(indices))
  pch_open <- rep(c(1, 2, 0),    length.out = length(indices))

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  # Bottom margin has to hold the axis, the axis label and the sub-caption
  # under it. mfrow shares one `mar` across both panels, so the two-panel case
  # takes the same allowance and its second panel simply has room to spare.
  if (isTRUE(entropy)) par(mfrow = c(2, 1), mar = c(5.4, 4.2, 2.6, 1))
  else                 par(mar = c(5.4, 4.2, 3, 1))

  Y <- as.matrix(tab[, indices, drop = FALSE])
  ylim <- range(Y[is.finite(Y)])
  plot(range(K), ylim, type = "n", xaxt = "n", xlab = "Number of classes",
       ylab = if (length(indices) == 1L) indices else "Information criterion",
       main = main %||% "Information criteria by number of classes")
  axis(1, at = K)

  for (j in seq_along(indices)) {
    y <- Y[, j]
    lines(K, y, col = cols[j], lty = ltys[j], lwd = 1.8)
    # An unreplicated maximum is not a fitted value to be trusted at face
    # value, so it is drawn hollow rather than solid.
    points(K, y, col = cols[j], pch = ifelse(unrep, pch_open[j], pch_fill[j]),
           cex = 1.2, lwd = 1.5)
    if (any(is.finite(y))) {
      i <- which.min(y)
      points(K[i], y[i], col = cols[j],
             pch = if (unrep[i]) pch_open[j] else pch_fill[j], cex = 2.1,
             lwd = 1.8)
      # Only the first index gets a guide line. One per index would turn the
      # panel into a grid and imply the criteria are voting.
      if (j == 1L)
        abline(v = K[i], lty = 3, col = adjustcolor(cols[j], alpha.f = 0.5))
    }
  }

  if (length(indices) > 1L)
    legend("topright", legend = indices, col = cols, lty = ltys,
           pch = pch_fill, lwd = 1.8, bty = "n", cex = 0.85)

  sub <- "Larger symbol marks the minimum; read the elbow, not the minimum alone."
  if (any(unrep))
    sub <- paste(sub, "Hollow symbols: maximum found by a single start.")
  mtext(sub, side = 1, line = 3.9, cex = 0.72, col = "grey30")

  if (isTRUE(entropy)) {
    ent <- tab$Entropy
    if (is.null(ent)) {
      warning("The fit table has no Entropy column; the second panel is empty.",
              call. = FALSE)
      ent <- rep(NA_real_, length(K))
    }
    # Fixed 0-1 limits. Autoscaling would magnify a range of a few hundredths
    # into a dramatic-looking curve.
    plot(range(K), c(0, 1), type = "n", xaxt = "n", ylim = c(0, 1),
         xlab = "Number of classes", ylab = "Relative entropy",
         main = "Classification certainty")
    axis(1, at = K)
    lines(K, ent, col = .okabe_ito[3], lty = 1, lwd = 1.8)
    points(K, ent, col = .okabe_ito[3],
           pch = ifelse(unrep, 1, 16), cex = 1.2, lwd = 1.5)
  }

  invisible(x)
}
