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

# Does a (possibly nested) measurement model contain any continuous block?
.has_continuous <- function(mm) {
  if (inherits(mm, "nested"))
    return(any(vapply(mm$models, .has_continuous, logical(1))))
  !is.null(mm$parameters$means)
}

# Build a [classes x items] matrix of values on a common [0, 1] scale.
# Column names carry a trailing "*" for any item that was rescaled. For nested
# (mixed-type) models the blocks are concatenated in submodel order.
.prepare_profile_data <- function(mm, indicators = NULL) {

  if (inherits(mm, "nested")) {
    blocks <- lapply(mm$models, .prepare_profile_data, indicators = indicators)
    return(do.call(cbind, blocks))
  }

  is_cont <- !is.null(mm$parameters$means)
  is_poly <- !is.null(mm$max_val)
  is_bin  <- !is.null(mm$parameters$pis) && is.null(mm$max_val)

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
  labels    <- if (is.null(class_labels)) paste("Class", seq_len(n_classes))
               else class_labels

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(6, 4, 4, 9), xpd = TRUE)

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

  legend(
    "topright",
    inset  = c(-0.24, 0),
    legend = labels,
    col    = my_colors,
    pch    = my_shapes,
    lty    = 1,
    lwd    = 2,
    bty    = "n",
    title  = "Class"
  )

  if (any(grepl("\\*", colnames(plot_mat))))
    mtext(
      "* Continuous items min-max scaled; ordinal items shown as scaled expected category.",
      side = 1, line = 4.5, adj = 0, cex = 0.8, font = 3, col = "grey30"
    )

  invisible(x)
}
