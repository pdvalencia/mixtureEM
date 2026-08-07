# ==============================================================================
# Repeated-Measures Latent Class Analysis (RMLCA / LLCA / LLPA)
# ==============================================================================

#' Repeated-Measures Latent Class Analysis
#'
#' @description
#' Fits a repeated-measures latent class model to data in which the same
#' indicators are observed at several occasions. Each person belongs to one
#' latent class for the whole study, so the classes describe *trajectories*:
#' patterns of response that span the occasions rather than a snapshot at one of
#' them. This is Collins and Lanza's RMLCA (2010, sec. 7.2); with continuous
#' indicators the same model is Wang and Wang's longitudinal latent profile
#' analysis (2020, sec. 6.3.1), and it is obtained here simply by setting
#' `measurement = "continuous"`.
#'
#' The likelihood is that of an ordinary latent class model applied to the
#' \eqn{J \times T} stacked indicators,
#' \deqn{P(y_i) = \sum_k \gamma_k \prod_t \prod_j \rho_{jtk}(y_{ijt}),}
#' so everything the package already offers — model selection, the bootstrap
#' likelihood-ratio test, predictors of class membership, distal outcomes,
#' survey designs and FIML for missing data — applies unchanged.
#'
#' `measurement_invariance` controls whether the item-response parameters
#' \eqn{\rho_{jtk}} are held equal across occasions. Constraining them makes a
#' class label mean the same thing at every occasion and sharply reduces the
#' number of free parameters; leaving them free lets an item behave differently
#' over time. The two models are nested, so [`lr_test()`] tests the
#' restriction directly.
#'
#' @param indicators The repeated indicators. Either a wide matrix or data frame
#'   with \eqn{J \times T} columns (see `layout`), a three-dimensional array with
#'   dimensions n by items by times, or a long data frame together with `id` and
#'   `time`.
#' @param n_classes Integer. Number of latent classes.
#' @param times Integer. Number of occasions. Required for wide input; inferred
#'   otherwise.
#' @param measurement Measurement model for a single occasion's items:
#'   `"binary"`, `"categorical"`, `"continuous"`, or a named list for a mixed
#'   block (as in [`fit_mixture()`]).
#' @param measurement_invariance Whether the item parameters are held equal
#'   across occasions. `"none"` (the default) estimates them separately at every
#'   occasion, which is usually what you want here: the classes are patterns of
#'   change, so forcing the items to behave identically over time can erase the
#'   very differences being modelled. `"full"` holds every item equal, and
#'   `"partial"` holds only the items named in `invariant_items`.
#' @param invariant_items Item indices or names held equal across occasions.
#'   Used only when `measurement_invariance = "partial"`.
#' @param layout For wide input, whether columns run `"time_major"`
#'   (all items of occasion 1, then all items of occasion 2, ...) or
#'   `"item_major"` (all occasions of item 1, then all occasions of item 2, ...).
#' @param id,time For long input, the case and occasion identifiers, given
#'   either as column names or as vectors.
#' @param items For long input, the columns to treat as indicators.
#' @param item_names,time_labels Optional display labels.
#' @param predictors Optional predictors of class membership. A grouping
#'   variable (Collins and Lanza, sec. 7.2.1) is entered this way: a
#'   multiple-group model and a model with the group as a dummy predictor are
#'   equivalent when measurement is invariant across groups (their sec. 6.10.2).
#' @param ... Further arguments passed to [`fit_mixture()`], such as `outcome`,
#'   `n_init`, `random_state`, `weights`, `strata` or `cluster`.
#'
#' @return An object of class `c("rmlca", "mixture_model")`. In addition to the
#'   usual fields it carries `$longitudinal`, holding the item and occasion
#'   labels, the invariance specification and the wave-missingness pattern.
#'
#' @references
#' Collins, L. M., & Lanza, S. T. (2010). \emph{Latent Class and Latent
#' Transition Analysis: With Applications in the Social, Behavioral, and Health
#' Sciences}. Wiley (sec. 6.10).
#'
#' @seealso [`fit_lta()`] for a model in which class membership may change
#'   between occasions, and [`lr_test()`] for testing invariance.
#' @export
fit_rmlca <- function(indicators,
                      n_classes = 2,
                      times = NULL,
                      measurement = "binary",
                      measurement_invariance = c("none", "full", "partial"),
                      invariant_items = NULL,
                      layout = c("time_major", "item_major"),
                      id = NULL, time = NULL, items = NULL,
                      item_names = NULL, time_labels = NULL,
                      predictors = NULL,
                      ...) {

  measurement_invariance <- match.arg(measurement_invariance)
  layout                 <- match.arg(layout)
  time_invariance        <- measurement_invariance

  prep <- .prepare_longitudinal(indicators, times = times, items = items,
                                layout = layout, id = id, time = time,
                                item_names = item_names,
                                time_labels = time_labels)

  spec <- .resolve_invariance(time_invariance, invariant_items,
                              prep$item_names, measurement)

  engine <- .longitudinal_measurement_spec(measurement, prep$X, prep$n_items,
                                           prep$n_times)

  if (!is.null(predictors))
    predictors <- .as_named_covariates(predictors, substitute(predictors),
                                       "predictor")

  fit <- fit_mixture(
    indicators      = prep$X,
    n_classes       = n_classes,
    measurement     = "time_blocks",
    predictors      = predictors,
    n_items         = prep$n_items,
    n_times         = prep$n_times,
    sub_model       = engine$sub_model,
    invariant_items = spec$invariant_items,
    max_val         = engine$max_val,
    ...
  )

  fit$longitudinal <- list(
    model           = "rmlca",
    n_items         = prep$n_items,
    n_times         = prep$n_times,
    item_names      = prep$item_names,
    time_labels     = prep$time_labels,
    time_invariance = time_invariance,
    invariant_items = spec$invariant_items,
    wave_missing    = prep$wave_missing,
    measurement     = measurement
  )
  class(fit) <- c("rmlca", class(fit))
  fit
}

# ------------------------------------------------------------------------------
# Shared helpers (also used by fit_lta)
# ------------------------------------------------------------------------------

# Turn the user's invariance request into an explicit vector of item indices.
.resolve_invariance <- function(time_invariance, invariant_items, item_names,
                                measurement) {
  J <- length(item_names)
  inv <- switch(
    time_invariance,
    none    = integer(0),
    full    = seq_len(J),
    partial = {
      if (is.null(invariant_items))
        stop('measurement_invariance = "partial" also needs `invariant_items`, ',
             "naming which items are held equal across occasions.",
             call. = FALSE)
      if (is.list(measurement))
        stop("Partial invariance is not defined for a mixed measurement block; ",
             'use "full" or "none".', call. = FALSE)
      idx <- if (is.character(invariant_items))
        match(invariant_items, item_names) else as.integer(invariant_items)
      if (anyNA(idx) || any(idx < 1L) || any(idx > J))
        stop("`invariant_items` must name or index items in the measurement ",
             "block.", call. = FALSE)
      sort(unique(idx))
    }
  )
  list(invariant_items = inv)
}

# Resolve the per-occasion measurement descriptor against the data: pick the
# FIML variant when anything is missing, and fix a single max_val for polytomous
# items so that every occasion shares the same category set (without which the
# invariance test would be comparing models over different response spaces).
.longitudinal_measurement_spec <- function(measurement, X, n_items, n_times) {
  sub_model <- .resolve_emission_descriptor(measurement, X)

  is_binary <- function(d) is.character(d) &&
    d %in% c("binary", "bernoulli", "binary_nan", "bernoulli_nan")
  is_poly   <- function(d) is.character(d) &&
    d %in% c("categorical", "multinoulli", "categorical_nan", "multinoulli_nan")

  if (is_binary(sub_model)) {
    vals <- X[!is.na(X)]
    if (length(vals) && !all(vals %in% c(0, 1)))
      stop('measurement = "binary" requires indicator values in {0, 1}.',
           call. = FALSE)
  }

  max_val <- NULL
  if (is_poly(sub_model)) {
    max_val <- max(X, na.rm = TRUE)
    if (!is.finite(max_val) || max_val != as.integer(max_val))
      stop('measurement = "categorical" requires integer-coded categories ',
           "(1, 2, 3, ...).", call. = FALSE)
    max_val <- as.integer(max_val)
  } else if (is.list(sub_model)) {
    # Mixed block: resolve each sub-block against the first occasion's columns
    # so that block-wise FIML upgrading matches the data actually seen there.
    sub_model <- .resolve_emission_descriptor(
      measurement, X[, .time_block_cols(1L, n_items), drop = FALSE])
  }

  list(sub_model = sub_model, max_val = max_val)
}

# Class-by-time-by-item array of the quantity that characterises each class at
# each occasion: endorsement probability (binary), expected category
# (polytomous), class mean (continuous) or event rate (count).
.rmlca_trajectories <- function(mm) {
  K  <- mm$n_components
  J  <- mm$n_items
  Tn <- mm$n_times
  arr  <- array(NA_real_, dim = c(K, Tn, J))
  kind <- "Probability"

  for (t in seq_len(Tn)) {
    sub <- mm$models[[t]]
    if (inherits(sub, "nested")) return(NULL)
    if (!is.null(sub$parameters$means)) {
      arr[, t, ] <- sub$parameters$means
      kind <- "Class mean"
    } else if (!is.null(sub$parameters$rates)) {
      arr[, t, ] <- sub$parameters$rates
      kind <- "Event rate"
    } else if (!is.null(sub$max_val)) {
      M <- sub$max_val
      for (j in seq_len(J)) {
        p <- sub$parameters$pis[, ((j - 1L) * M + 1L):(j * M), drop = FALSE]
        arr[, t, j] <- rowSums(sweep(p, 2, seq_len(M), "*"))
      }
      kind <- "Expected category"
    } else if (!is.null(sub$parameters$pis)) {
      arr[, t, ] <- sub$parameters$pis
    } else {
      return(NULL)
    }
  }
  list(values = arr, kind = kind)
}

# ------------------------------------------------------------------------------
# Methods
# ------------------------------------------------------------------------------

#' Print a Fitted Repeated-Measures Latent Class Model
#'
#' @param x An object returned by [`fit_rmlca()`].
#' @param ... Passed to the next method.
#' @return `x`, invisibly.
#' @export
print.rmlca <- function(x, ...) {
  lg <- x$longitudinal
  cat("\n")
  cat("=========================================================\n")
  cat("        REPEATED-MEASURES LATENT CLASS MODEL\n")
  cat("=========================================================\n")
  cat(sprintf("Items x Occasions  : %d x %d\n", lg$n_items, lg$n_times))
  cat(sprintf("Item parameters    : %s, %s\n",
              if (is.list(lg$measurement)) "mixed" else lg$measurement,
              switch(lg$time_invariance,
                     none    = "estimated separately at each occasion",
                     full    = "held equal across occasions",
                     partial = sprintf("%s held equal across occasions",
                                       paste(lg$item_names[lg$invariant_items],
                                             collapse = ", ")))))
  if (any(lg$wave_missing))
    cat(sprintf("Wave attrition     : %d case-occasions with no observed item\n",
                sum(lg$wave_missing)))
  NextMethod()
}

#' Trajectory Plot for a Repeated-Measures Latent Class Model
#'
#' @description
#' Draws one panel per item with occasions on the x-axis and one line per latent
#' class, which is the natural way to read RMLCA output: the classes *are* the
#' trajectories. Pass `type = "profile"` for the single-panel indicator profile
#' used by [`plot.mixture_model()`].
#'
#' @param x An object returned by [`fit_rmlca()`].
#' @param type `"trajectory"` (default) or `"profile"`.
#' @param main Plot title.
#' @param class_labels Optional class labels for the legend.
#' @param colors Optional colour vector, recycled across classes.
#' @param ... Passed to the profile plot when `type = "profile"`.
#' @return `x`, invisibly.
#' @importFrom graphics par matplot axis legend mtext plot.new
#' @export
plot.rmlca <- function(x, type = c("trajectory", "profile"),
                       main = NULL, class_labels = NULL, colors = NULL, ...) {
  type <- match.arg(type)
  if (type == "profile") {
    return(plot.mixture_model(
      x, main = main %||% "Latent Class / Profile Plot",
      class_labels = class_labels, colors = colors, ...))
  }

  traj <- .rmlca_trajectories(x$mm)
  if (is.null(traj))
    stop("Trajectories are not defined for this measurement model; ",
         'use type = "profile".', call. = FALSE)

  lg  <- x$longitudinal
  K   <- x$n_components
  Tn  <- lg$n_times
  J   <- lg$n_items
  cols   <- if (is.null(colors)) rep(.okabe_ito, length.out = K)
            else rep(colors, length.out = K)
  shapes <- rep(15:20, length.out = K)
  base   <- if (is.null(class_labels)) paste("Class", seq_len(K)) else class_labels
  labels <- .class_plot_labels(base, x$weights)

  ylim <- if (traj$kind == "Class mean") range(traj$values, na.rm = TRUE)
          else if (traj$kind == "Expected category")
            range(c(1, traj$values), na.rm = TRUE)
          else c(0, 1)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  nr <- ceiling(sqrt(J)); nc <- ceiling(J / nr)

  # The legend never goes inside a panel, where it would cover data. When the
  # grid leaves an empty cell it goes there, which is both readable and free;
  # otherwise it gets a reserved strip in the bottom outer margin, deep enough
  # for however many rows the class labels need.
  spare    <- nr * nc > J
  leg_cols <- if (spare) 1L else max(1L, min(K, 3L))
  leg_rows <- ceiling(K / leg_cols)
  par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1),
      oma = c(if (spare) 0 else 2.5 + 1.1 * leg_rows, 0, 3, 0))

  for (j in seq_len(J)) {
    vals_j <- matrix(traj$values[, , j], nrow = K, ncol = Tn)
    matplot(seq_len(Tn), t(vals_j),
            type = "b", pch = shapes, lty = 1, lwd = 2, col = cols,
            ylim = ylim, xaxt = "n", xlab = "", ylab = traj$kind,
            main = lg$item_names[j], bty = "l", las = 1)
    axis(1, at = seq_len(Tn), labels = lg$time_labels)
  }

  mtext(main %||% "Latent class trajectories", outer = TRUE, line = 1,
        cex = 1.1, font = 2)

  if (spare) {
    plot.new()
    legend("center", legend = labels, col = cols, pch = shapes, lty = 1,
           lwd = 2, bty = "n", ncol = leg_cols, cex = 1)
  } else {
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
        new = TRUE)
    plot.new()
    legend("bottom", legend = labels, col = cols, pch = shapes, lty = 1,
           lwd = 2, bty = "n", ncol = leg_cols, cex = 0.85, xpd = TRUE)
  }

  invisible(x)
}
