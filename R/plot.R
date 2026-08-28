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

# A one-line footnote for a fit whose class prevalence varies by group
# (`group_effects = "prevalence"` or `"both"`): the legend percentages above
# come from `x$weights`, the pooled marginal shares across all groups, not the
# per-group breakdown the model actually estimates. NULL for every other fit,
# so the mtext() call this feeds is silent for the common case.
.group_pooled_note <- function(x) {
  if (is.null(.group_gamma_matrix(x))) return(NULL)
  "Legend percentages are pooled across groups; see class_sizes() for the by-group breakdown."
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

# ------------------------------------------------------------------------------
# Standardized profile heights for the bar chart
# ------------------------------------------------------------------------------
# A [classes x items] matrix of z-scored conditional means. The profile plot's
# min-max scaling is hostage to a single extreme observation and has no
# meaningful origin; a z-score has both a scale-free unit and a zero that means
# "the average case". The property that earns the extra code path is that the
# picture is unchanged by standardizing the indicators before fitting, so
# standardization becomes a display choice rather than a data-preparation step.
#
# Denominators:
#   "total"  -> the weighted observed SD of the column, which is what "we
#               standardized the indicators" means in the applied literature
#               and is what makes the figure invariant to pre-standardization.
#   "within" -> the model-implied within-class SD, pooled over classes by their
#               weights. Reads as a Cohen's-d-like effect against residual
#               rather than total dispersion, and is the more informative of the
#               two; it is not the default only because it is not what applied
#               papers plot.
#
# `type` names the renderer that asked, so the two error messages below name the
# argument the user actually passed rather than always naming the bar chart.
.profile_bar_heights <- function(x, scale = "total", type = "bar") {

  mm <- x$mm
  if (!(class(mm)[1] %in% c("gaussian_diag", "gaussian_diag_nan",
                            "gaussian_unit", "gaussian_unit_nan")))
    stop("plot(type = \"", type, "\") standardizes conditional means, so it ",
         "needs an all-continuous measurement model. This fit has ",
         class(mm)[1],
         " emissions; use type = \"profile\", which places indicators of any ",
         "type on a common axis.", call. = FALSE)

  if (is.null(x$data))
    stop("Raw indicators are not stored on this fit, so conditional means ",
         "cannot be standardized against the observed data. Refit with the ",
         "current version, or use type = \"profile\".", call. = FALSE)
  # (`type` is unused past this point; the heights themselves do not depend on
  # which renderer will draw them, which is the reason there is one helper.)

  means <- mm$parameters$means
  cols  <- colnames(means) %||% paste0("Cont_", seq_len(ncol(means)))
  # Match the observed columns by name, exactly as .prepare_profile_data() does.
  # An unnamed indicator matrix leaves the means unnamed too, and there the
  # columns are in fitting order and can only be matched positionally.
  matched <- .match_indicator_columns(cols, x$data)
  if (is.null(matched))
    stop("Could not match the class means to the stored indicators by name, ",
         "so they cannot be standardized. Use type = \"profile\".",
         call. = FALSE)

  w  <- x$sample_weights %||% rep(1, nrow(x$data))
  H  <- means
  for (j in seq_len(ncol(means))) {
    obs <- x$data[, matched[j]]
    ok  <- !is.na(obs)
    wj  <- w[ok]; oj <- obs[ok]
    m_j <- sum(wj * oj) / sum(wj)

    s_j <- if (identical(scale, "within")) {
      # Unit-variance models fix the within-class variance at 1 by definition.
      if (is.null(mm$parameters$covariances)) 1
      else sqrt(sum(x$weights * mm$parameters$covariances[, j]))
    } else {
      sqrt(sum(wj * (oj - m_j)^2) / sum(wj))
    }
    if (!is.finite(s_j) || s_j <= 0) s_j <- 1  # a constant column: leave centred

    H[, j] <- (means[, j] - m_j) / s_j
  }

  colnames(H) <- cols
  H
}

#' Profile Plot for a Fitted Mixture Model
#'
#' @description
#' Draws the measurement model in one of three ways.
#'
#' `type = "profile"` (the default) is a line plot: one line per latent class,
#' with every indicator placed on a common \[0, 1\] axis. Binary indicators are
#' shown as endorsement probabilities; continuous indicators are min-max scaled
#' against their observed range; polytomous indicators are summarised by their
#' expected category and scaled to \[0, 1\]. Rescaled items are marked with "*".
#'
#' `type = "bar"` is the grouped bar chart of standardized class means that
#' applied latent-profile papers publish, and requires an all-continuous
#' measurement model. Indicators sit on the x-axis with one bar per class inside
#' each group; a bar above the zero line is an above-average conditional mean
#' and one below it a below-average mean. Heights are z-scores rather than the
#' profile plot's min-max scaling, which is hostage to a single extreme
#' observation and has no meaningful origin. Because a z-score is scale-free,
#' the figure is the same whether or not the indicators were standardized before
#' fitting, so standardization is a display choice here rather than a step in
#' preparing the data.
#'
#' `type = "line"` draws those same z-scores as one connected line per class,
#' with a zero reference line, and likewise requires an all-continuous
#' measurement model. It answers a different question from the bar chart off the
#' same numbers: bars group by indicator and invite comparing classes one
#' indicator at a time, while a line follows a single class across all of them
#' and shows the *shape* of its profile — whether two classes differ in level or
#' in pattern. Prefer it over `type = "profile"` whenever every indicator is
#' continuous, since the min-max axis there has no meaningful origin and its
#' shape depends on the sample's most extreme observation.
#'
#' Uses only base graphics and the colour-blind-friendly Okabe-Ito palette.
#'
#' @param x A fitted \code{mixture_model} object.
#' @param type One of `"profile"` (the default min-max line plot), `"bar"` (the
#'   standardized profile bar chart) or `"line"` (the same standardized profile
#'   drawn as lines), all described above. `"bar"` and `"line"` require an
#'   all-continuous measurement model.
#' @param main Plot title. Defaults to a title chosen for `type`.
#' @param class_labels Optional character vector of labels for the classes.
#'   Defaults to "Class 1", "Class 2", ...
#' @param colors Optional vector of colours (one per class). Defaults to the
#'   Okabe-Ito palette, recycled if necessary.
#' @param scale For `type = "bar"` and `type = "line"`, the denominator of the
#'   z-score.
#'   `"total"` (the default) divides by the weighted observed SD of the
#'   indicator, which is what standardizing the indicators means in the applied
#'   literature and is what makes the figure invariant to pre-standardization.
#'   `"within"` divides by the model-implied within-class SD pooled over
#'   classes, giving a Cohen's-d-like reading against residual rather than total
#'   dispersion. `"within"` is arguably the more informative quantity; it is not
#'   the default because it is not what applied papers plot.
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
#' @importFrom graphics par matplot axis text legend mtext barplot abline
#' @export
plot.mixture_model <- function(x, type = c("profile", "bar", "line"),
                               main = NULL,
                               class_labels = NULL, colors = NULL,
                               scale = c("total", "within"), ...) {

  if (is.null(x$mm)) stop("No measurement model to plot.")

  type  <- match.arg(type)
  scale <- match.arg(scale)

  if (identical(type, "bar"))
    return(.plot_profile_bar(x, main = main, class_labels = class_labels,
                             colors = colors, scale = scale))

  if (identical(type, "line"))
    return(.plot_profile_line(x, main = main, class_labels = class_labels,
                              colors = colors, scale = scale))

  if (is.null(main)) main <- "Latent Class / Profile Plot"

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
  n_footnotes <- as.integer(any(grepl("\\*", colnames(plot_mat)))) +
    as.integer(!is.null(.group_pooled_note(x)))
  mar_bottom <- max(6, min(3.5 + 0.4 * max(nchar(colnames(plot_mat))), 12)) +
    max(0, n_footnotes - 1)
  # The first item's rotated label leans left of the axis, so the left margin
  # grows with it too (the y-axis title needs the 4-line floor regardless).
  mar_left   <- max(4, min(2 + 0.25 * nchar(colnames(plot_mat)[1]), 8))
  # A `main` with embedded line breaks needs the top margin to grow with it --
  # base R's title() does not reserve extra room for a multi-line string on
  # its own.
  n_title_lines <- lengths(gregexpr("\n", main %||% "")) + 1
  mar_top <- 2 + 2 * n_title_lines
  par(mar = c(mar_bottom, mar_left, mar_top, mar_right), xpd = TRUE)

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

  footnotes <- character(0)
  if (any(grepl("\\*", colnames(plot_mat))))
    footnotes <- c(footnotes,
                   "* Continuous items min-max scaled; ordinal items shown as scaled expected category.")
  footnotes <- c(footnotes, .group_pooled_note(x))
  for (i in seq_along(footnotes))
    mtext(
      footnotes[i],
      side = 1, line = mar_bottom - 1.2 + (i - 1), adj = 0, cex = 0.8, font = 3,
      col = "grey30"
    )

  invisible(x)
}

# The grouped bar chart behind plot(type = "bar"). Groups are indicators and
# bars within a group are classes, which is the arrangement applied latent
# profile papers use: the reader compares classes on one indicator at a time.
# barplot(beside = TRUE) reads a matrix as one group per *column* and one bar
# per *row*, so the [classes x items] matrix goes in untransposed. It draws a
# negative height downward from zero on its own, so the only things to add are
# the zero line itself and a symmetric ylim.
.plot_profile_bar <- function(x, main = NULL, class_labels = NULL,
                              colors = NULL, scale = "total") {

  H <- .profile_bar_heights(x, scale = scale)

  n_classes <- nrow(H)
  n_items   <- ncol(H)

  if (is.null(main))
    main <- "Standardized Class Profiles"

  my_colors <- if (is.null(colors)) rep(.okabe_ito, length.out = n_classes)
               else rep(colors, length.out = n_classes)

  base_labels <- if (is.null(class_labels)) paste("Class", seq_len(n_classes))
                 else class_labels
  labels <- .class_plot_labels(base_labels, x$weights)

  item_labels <- .shorten_labels(colnames(H), width = 16L)
  # Rotate only when the names are long enough to collide; indicators are often
  # named "a" or "x1", and a 45-degree one-character label just looks broken.
  rotate <- max(nchar(item_labels)) > 8L

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  # Same reasoning as the profile plot: the right margin is sized by the legend
  # it holds and the bottom margin by the item names, which need room only in
  # proportion to their length when they are rotated.
  gp_note <- .group_pooled_note(x)
  mar_right  <- min(3 + 0.45 * max(nchar(labels)), 18)
  mar_bottom <- (if (rotate)
    max(6, min(3.5 + 0.4 * max(nchar(item_labels)), 12)) else 4) +
    (if (is.null(gp_note)) 0 else 1)
  par(mar = c(mar_bottom, 4.5, 4, mar_right), xpd = TRUE)

  lim  <- max(abs(H), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  ylim <- c(-1.1 * lim, 1.1 * lim)

  y_lab <- if (identical(scale, "within"))
    "Class mean (SD within class)" else "Class mean (SD)"

  at <- barplot(
    H, beside = TRUE, col = my_colors, border = NA, names.arg = rep("", n_items),
    ylim = ylim, ylab = y_lab, main = main, las = 1
  )

  abline(h = 0, col = "grey20")

  # One label per indicator group, centred under its cluster of bars. The offset
  # is measured from the plotting region barplot actually established, not from
  # the requested ylim -- barplot extends the region past it, so an ylim-based
  # offset drops the labels well below the axis.
  centres <- colMeans(at)
  text(
    x      = centres,
    y      = par("usr")[3] - 0.03 * diff(par("usr")[3:4]),
    labels = item_labels,
    srt    = if (rotate) 45 else 0,
    adj    = if (rotate) 1 else c(0.5, 1),
    cex    = 0.9
  )

  x_pos <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  legend(
    x = x_pos, y = par("usr")[4], legend = labels, fill = my_colors,
    border = NA, bty = "n", cex = 0.9, title = "Class"
  )

  if (!is.null(gp_note))
    mtext(gp_note, side = 1, line = mar_bottom - 1.2, adj = 0, cex = 0.8,
          font = 3, col = "grey30")

  invisible(x)
}

# The line plot behind plot(type = "line"): the bar chart's z-scores, drawn as
# one connected line per class instead of one bar per class.
#
# It exists alongside the bar chart because the two answer different questions
# off the same numbers. Bars group by indicator and invite comparing classes
# within an indicator; a line follows one class across all of them and shows the
# *shape* of its profile - which is what "profile" means in the applied
# literature, and what a reader looks for when asking whether two classes differ
# in level or in pattern. The default type = "profile" already draws lines, but
# on a min-max axis, where the shape is hostage to a single extreme observation
# and the vertical position has no meaning beyond rank within the sample range.
#
# Two things differ from the min-max renderer. The y-axis is symmetric about
# zero rather than [0, 1], and zero is drawn, because on a z-score scale it is
# the sample average and every reading is relative to it. There are no "*"
# suffixes on the item names: every column is standardized the same way, so that
# is a fact about the axis, stated once in its label, not a mark on each item.
.plot_profile_line <- function(x, main = NULL, class_labels = NULL,
                               colors = NULL, scale = "total") {

  H <- .profile_bar_heights(x, scale = scale, type = "line")

  n_classes <- nrow(H)
  n_items   <- ncol(H)

  if (is.null(main)) main <- "Standardized Class Profiles"

  my_colors <- if (is.null(colors)) rep(.okabe_ito, length.out = n_classes)
               else rep(colors, length.out = n_classes)
  my_shapes <- rep(15:20, length.out = n_classes)

  base_labels <- if (is.null(class_labels)) paste("Class", seq_len(n_classes))
                 else class_labels
  labels <- .class_plot_labels(base_labels, x$weights)

  item_labels <- .shorten_labels(colnames(H), width = 30L)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  gp_note <- .group_pooled_note(x)
  mar_right  <- min(3 + 0.45 * max(nchar(labels)), 18)
  mar_bottom <- max(6, min(3.5 + 0.4 * max(nchar(item_labels)), 12)) +
    (if (is.null(gp_note)) 0 else 1)
  mar_left   <- max(4.5, min(2 + 0.25 * nchar(item_labels[1]), 8))
  par(mar = c(mar_bottom, mar_left, 4, mar_right), xpd = TRUE)

  lim <- max(abs(H), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  y_lab <- if (identical(scale, "within"))
    "Class mean (SD within class)" else "Class mean (SD)"

  matplot(
    x    = seq_len(n_items),
    y    = t(H),
    type = "b",
    pch  = my_shapes,
    lty  = 1,
    lwd  = 2,
    col  = my_colors,
    ylim = c(-1.1 * lim, 1.1 * lim),
    xlim = c(0.75, n_items + 0.25),
    xaxt = "n",
    xlab = "",
    ylab = y_lab,
    main = main,
    bty  = "l",
    las  = 1
  )

  abline(h = 0, col = "grey20")

  axis(1, at = seq_len(n_items), labels = FALSE)
  text(
    x      = seq_len(n_items),
    y      = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
    labels = item_labels,
    srt    = 45,
    adj    = 1,
    cex    = 0.9
  )

  x_pos <- par("usr")[2] + 0.02 * diff(par("usr")[1:2])
  legend(
    x = x_pos, y = par("usr")[4], legend = labels, col = my_colors,
    pch = my_shapes, lty = 1, lwd = 2, bty = "n", cex = 0.9, title = "Class"
  )

  if (!is.null(gp_note))
    mtext(gp_note, side = 1, line = mar_bottom - 1.2, adj = 0, cex = 0.8,
          font = 3, col = "grey30")

  .cat_label_legend(item_labels, indent = "")

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
