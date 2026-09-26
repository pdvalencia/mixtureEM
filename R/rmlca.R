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
#'   `n_init`, `random_state`, `weights`, `strata` or `cluster`. The EM
#'   convergence rule is [`fit_mixture()`]'s fixed one, since this function
#'   estimates through it; there is no `tol`-style argument to pass here.
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
  prep$X <- engine$X
  .warn_empty_categories(engine$empty, prep$item_names, spec$invariant_items,
                         prep$time_labels, "occasion")

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
    cats            = engine$cats,
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
#
# `measurement` is accepted (rather than dropped from the signature) because a
# mixed measurement block changes what an *integer* `invariant_items` means:
# .normalize_measurement() regroups indicator columns by type, so item_names
# (and therefore a name lookup below) already reflects the regrouped order,
# but a caller naming items by their original integer position does not get
# that translation. A caller documents this in its own @param; nothing here
# needs to special-case `measurement` beyond accepting it for that note.
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

# Each categorical item's number of categories, pooled over every block (group
# or occasion) of a J-item layout, and the cells with no response in them.
#
# The categories an item can take are a property of the item, not of the
# sample one group or one occasion happens to supply. Counting them per block
# would give the configural and the invariant model different response spaces,
# and their likelihoods could no longer be compared. So an item's count is its
# highest code anywhere, or `declared` -- a factor's full set of levels -- when
# that is higher, and a category a block never uses is estimated at 0 there
# rather than removed. `items` are positions in 1..n_items; `declared` is
# aligned with `items`. Returned: `cats` (aligned with `items`), `pooled` (item,
# category) with no response in any block, and `by_block` (item, block,
# category) with none in that block although the item was observed there and
# the category is used elsewhere.
.pooled_block_cats <- function(X, items, n_items, n_blocks, declared = NULL) {
  cats     <- integer(length(items))
  pooled   <- list()
  by_block <- list()
  for (i in seq_along(items)) {
    j    <- items[i]
    cols <- (seq_len(n_blocks) - 1L) * n_items + j
    v    <- X[, cols]
    v    <- v[!is.na(v)]
    d    <- if (is.null(declared) || is.na(declared[i])) 0L else as.integer(declared[i])
    cats[i] <- max(2L, d, if (length(v)) as.integer(max(v)) else 0L)
    # Codes that are not 1, 2, 3, ... are refused with their own message when
    # the emission is fitted; reporting empty categories for them first would
    # only bury it.
    if (!length(v) || any(v < 1) || any(v != floor(v))) next
    zero <- which(tabulate(v, nbins = cats[i]) == 0L)
    if (length(zero))
      pooled[[length(pooled) + 1L]] <- data.frame(item = j, category = zero)
    if (n_blocks < 2L) next
    for (b in seq_len(n_blocks)) {
      vb <- X[, cols[b]]
      vb <- vb[!is.na(vb)]
      if (!length(vb)) next
      zb <- setdiff(which(tabulate(vb, nbins = cats[i]) == 0L), zero)
      if (length(zb))
        by_block[[length(by_block) + 1L]] <-
          data.frame(item = j, block = b, category = zb)
    }
  }
  list(cats     = cats,
       pooled   = if (length(pooled)) do.call(rbind, pooled) else NULL,
       by_block = if (length(by_block)) do.call(rbind, by_block) else NULL)
}

# Say which categorical items have a category nobody used.
#
# Neither case is an error. A category with no response in the whole sample --
# an unused factor level, or a gap in the codes -- is estimated at 0 in every
# class and still counted, because it is part of the item's declared response
# space. A category with no response in one group or occasion only matters
# where that block's probabilities are free: there it is a boundary estimate
# with no standard error. An item held equal across the blocks pools every
# block's responses, so it is left out of that half of the report.
.warn_empty_categories <- function(empty, item_names, invariant_items = integer(0),
                                   block_labels = NULL, block_word = "group",
                                   max_show = 6L) {
  if (is.null(empty)) return(invisible(NULL))
  nm <- function(j) item_names[j] %||% paste0("item ", j)
  show <- function(lines) {
    more <- length(lines) - max_show
    paste0(paste(utils::head(lines, max_show), collapse = "; "),
           if (more > 0L) sprintf("; and %d more", more) else "")
  }

  p <- empty$pooled
  if (!is.null(p) && nrow(p)) {
    lines <- vapply(split(p, p$item), function(d)
      sprintf("%s (category %s)", nm(d$item[1L]), paste(d$category, collapse = ", ")),
      character(1))
    warning(paste0(
      "Some categorical items have a category with no responses: ", show(lines),
      ". Each is estimated at 0 in every class and still counted as a ",
      "parameter. If the category is not part of the item, recode it (or drop ",
      "the unused factor level) and refit."), call. = FALSE)
  }

  b <- empty$by_block
  if (!is.null(b) && nrow(b)) b <- b[!b$item %in% invariant_items, , drop = FALSE]
  if (!is.null(b) && nrow(b)) {
    lab <- function(k) block_labels[k] %||% as.character(k)
    lines <- vapply(split(b, list(b$item, b$block), drop = TRUE), function(d)
      sprintf("%s in %s %s (category %s)", nm(d$item[1L]), block_word,
              lab(d$block[1L]), paste(d$category, collapse = ", ")),
      character(1))
    warning(paste0(
      "Some categories have no responses in one ", block_word, ": ", show(lines),
      ". Where the item's probabilities differ by ", block_word, ", each of ",
      "these is estimated at 0, a boundary estimate with no standard error. ",
      "Holding the item equal across ", block_word, "s pools the responses."),
      call. = FALSE)
  }
  invisible(NULL)
}

# Resolve the per-occasion measurement descriptor against the data: pick the
# FIML variant when anything is missing, and fix each categorical item's number
# of categories over every occasion, so all occasions share one response space
# (without which the invariance test would be comparing models over different
# response spaces). `declared` optionally gives each item's declared count, a
# factor's number of levels, as a floor. `empty` reports the categories nobody
# used, for .warn_empty_categories().
.longitudinal_measurement_spec <- function(measurement, X, n_items, n_times,
                                           declared = NULL) {
  sub_model <- .resolve_emission_descriptor(measurement, X)

  is_binary  <- function(d) is.character(d) &&
    d %in% c("binary", "bernoulli", "binary_nan", "bernoulli_nan")
  is_poly    <- function(d) is.character(d) &&
    d %in% c("categorical", "multinoulli", "categorical_nan", "multinoulli_nan")
  is_ordinal <- function(d) is.character(d) && d %in% c("ordinal", "ordinal_nan")

  if (is_binary(sub_model)) {
    X <- .recode_binary_blocks(X, n_items, n_times)$X
    vals <- X[!is.na(X)]
    if (length(vals) && !all(vals %in% c(0, 1)))
      stop('measurement = "binary" requires indicator values in {0, 1}, or a ',
           "two-level coding of each item that can be mapped to them.",
           call. = FALSE)
  }

  max_val <- NULL
  cats    <- NULL
  empty   <- NULL
  if (is_poly(sub_model)) {
    max_val <- max(X, na.rm = TRUE)
    if (!is.finite(max_val) || max_val != as.integer(max_val))
      stop('measurement = "categorical" requires integer-coded categories ',
           "(1, 2, 3, ...).", call. = FALSE)
    empty   <- .pooled_block_cats(X, seq_len(n_items), n_items, n_times, declared)
    cats    <- empty$cats
    max_val <- max(cats)
  } else if (is_ordinal(sub_model)) {
    # Per-item category counts, read across every occasion so all T columns
    # of one item share the same response space (### 14.18.6). Codes must be
    # 1-based and contiguous, and no interior category may be empty: a
    # category with zero observed responses that is neither the top nor the
    # bottom one leaves one of its thresholds unidentified.
    cats <- integer(n_items)
    for (j in seq_len(n_items)) {
      cols <- ((seq_len(n_times) - 1L) * n_items) + j
      vals <- X[, cols][!is.na(X[, cols])]
      if (!length(vals)) {
        cats[j] <- 2L
        next
      }
      if (!all(vals == as.integer(vals)) || min(vals) < 1)
        stop(sprintf(paste0(
          'measurement = "ordinal" requires integer-coded categories ',
          "(1, 2, 3, ...) for item %d."), j), call. = FALSE)
      cats[j] <- as.integer(max(vals))
      observed <- tabulate(vals, nbins = cats[j])
      interior <- if (cats[j] > 2L) seq(2L, cats[j] - 1L) else integer(0)
      if (any(observed[interior] == 0L))
        stop(sprintf(paste(
          "measurement = \"ordinal\": item %d has an empty interior",
          "category (a category between the lowest and highest observed",
          "value with zero responses), which leaves one of its thresholds",
          "unidentified."), j), call. = FALSE)
    }
  } else if (is.list(sub_model)) {
    # Mixed block: resolve each sub-block against the first occasion's columns
    # so that block-wise FIML upgrading matches the data actually seen there.
    sub_model <- .resolve_emission_descriptor(
      measurement, X[, .time_block_cols(1L, n_items), drop = FALSE])
    # A categorical sub-block's items get their counts pooled over occasions
    # here too; left to each occasion's own data, a category one occasion
    # happens not to use would change that occasion's response space.
    offset <- 0L
    for (key in names(sub_model)) {
      n_j <- sub_model[[key]]$n_columns
      if (is_poly(sub_model[[key]]$model)) {
        items <- offset + seq_len(n_j)
        pc <- .pooled_block_cats(X, items, n_items, n_times,
                                 sub_model[[key]]$cats %||% declared[items])
        sub_model[[key]]$cats <- pc$cats
        empty <- list(pooled   = rbind(empty$pooled, pc$pooled),
                      by_block = rbind(empty$by_block, pc$by_block))
      }
      offset <- offset + n_j
    }
  }

  list(sub_model = sub_model, max_val = max_val, cats = cats, X = X,
       empty = empty)
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
