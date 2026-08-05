# ==============================================================================
# Longitudinal data preparation
# ==============================================================================
#
# Both longitudinal models in this package see the same internal layout: an
# n x (J * T) matrix in TIME-MAJOR order, i.e. all J items of time 1, then all J
# items of time 2, and so on. Repeated-measures LCA consumes it as one long row
# of indicators partitioned into per-time blocks; latent transition analysis
# slices one time block at a time and hands each to an ordinary measurement
# model. Everything downstream assumes this order, so all reshaping happens here.

# Build the canonical "item@time" column labels used throughout the longitudinal
# summaries and plots.
.longitudinal_colnames <- function(item_names, time_labels) {
  as.vector(vapply(time_labels,
                   function(tt) paste0(item_names, "@", tt),
                   character(length(item_names))))
}

# Validate and reshape user input into the canonical time-major matrix.
#
# `indicators` may be
#   * a wide matrix / data frame with J * T columns (see `layout`),
#   * a 3-dimensional array with dimensions n x J x T, or
#   * a long data frame, in which case `id` and `time` must be supplied and
#     every remaining selected column is treated as an item.
#
# `layout` describes how the columns of a wide matrix are ordered:
#   "time_major" - t1i1, t1i2, ..., t2i1, t2i2, ...  (blocks of items per time)
#   "item_major" - i1t1, i1t2, ..., i2t1, i2t2, ...  (blocks of times per item)
#
# Returns the canonical matrix plus the metadata every caller needs: item and
# time labels, and a logical n x T matrix flagging waves where a case has no
# observed item at all (complete wave attrition), which is what makes the FIML
# path explicit rather than incidental.
.prepare_longitudinal <- function(indicators, times = NULL, items = NULL,
                                  layout = c("time_major", "item_major"),
                                  id = NULL, time = NULL,
                                  item_names = NULL, time_labels = NULL) {
  layout <- match.arg(layout)

  # ---- 3-d array input -------------------------------------------------------
  if (!is.null(dim(indicators)) && length(dim(indicators)) == 3L) {
    d <- dim(indicators)
    n <- d[1]; J <- d[2]; Tn <- d[3]
    X <- matrix(indicators, nrow = n, ncol = J * Tn)  # column-major: j fastest
    if (is.null(item_names))
      item_names <- dimnames(indicators)[[2]] %||% paste0("Item", seq_len(J))
    if (is.null(time_labels))
      time_labels <- dimnames(indicators)[[3]] %||% paste0("T", seq_len(Tn))

  # ---- long data frame input -------------------------------------------------
  } else if (!is.null(id) || !is.null(time)) {
    if (is.null(id) || is.null(time))
      stop("Long-format input requires both `id` and `time`.", call. = FALSE)
    df <- as.data.frame(indicators)
    reserved <- character(0)
    if (length(id) == 1L && is.character(id)) {
      reserved <- c(reserved, id);   id   <- df[[id]]
    }
    if (length(time) == 1L && is.character(time)) {
      reserved <- c(reserved, time); time <- df[[time]]
    }
    if (length(id) != nrow(df) || length(time) != nrow(df))
      stop("`id` and `time` must have one entry per row of `indicators`.",
           call. = FALSE)

    item_cols <- if (!is.null(items)) items else setdiff(names(df), reserved)
    if (is.character(item_cols)) {
      miss <- setdiff(item_cols, names(df))
      if (length(miss))
        stop("Unknown item column(s): ", paste(miss, collapse = ", "),
             call. = FALSE)
    }
    item_mat <- data.matrix(df[, item_cols, drop = FALSE])

    ids   <- unique(id)
    tvals <- sort(unique(time))
    n <- length(ids); J <- ncol(item_mat); Tn <- length(tvals)

    X <- matrix(NA_real_, nrow = n, ncol = J * Tn)
    row_of  <- match(id, ids)
    time_of <- match(time, tvals)
    for (t in seq_len(Tn)) {
      sel <- which(time_of == t)
      if (length(sel))
        X[row_of[sel], ((t - 1) * J + 1):(t * J)] <- item_mat[sel, , drop = FALSE]
    }
    if (is.null(item_names))  item_names  <- colnames(item_mat)
    if (is.null(time_labels)) time_labels <- as.character(tvals)

  # ---- wide matrix input -----------------------------------------------------
  } else {
    X <- if (is.data.frame(indicators)) data.matrix(indicators) else
      as.matrix(indicators)
    if (is.null(times))
      stop("`times` must be given for wide-format input.", call. = FALSE)
    Tn <- as.integer(times)
    if (is.na(Tn) || Tn < 2L)
      stop("`times` must be an integer >= 2.", call. = FALSE)
    if (ncol(X) %% Tn != 0L)
      stop(sprintf(
        "The number of indicator columns (%d) is not a multiple of times (%d).",
        ncol(X), Tn), call. = FALSE)
    J <- ncol(X) %/% Tn

    orig_names <- colnames(X)
    if (layout == "item_major") {
      # i1t1, i1t2, ..., i2t1, ... -> reorder into time-major blocks
      idx <- as.vector(sapply(seq_len(Tn),
                              function(t) seq(t, by = Tn, length.out = J)))
      X <- X[, idx, drop = FALSE]
      if (!is.null(orig_names)) orig_names <- orig_names[idx]
    }
    if (is.null(item_names)) {
      raw <- if (!is.null(orig_names)) orig_names[seq_len(J)] else
        paste0("Item", seq_len(J))
      # Drop the occasion marker the columns carry, so that "drug_t1" and
      # "drug_t2" both display as "drug". A separator is required, otherwise
      # "item1" would be stripped down to nothing; the marker itself may be any
      # short label ("_t1", "_v1", ".wave1") or absent ("_1").
      item_names <- sub("[._ -][A-Za-z]{0,4}0*1$", "", raw)
      if (any(duplicated(item_names)) || any(!nzchar(item_names)))
        item_names <- raw
    }
    if (is.null(time_labels)) time_labels <- paste0("T", seq_len(Tn))
  }

  if (length(item_names) != J)
    stop("`item_names` must have one entry per item.", call. = FALSE)
  if (length(time_labels) != Tn)
    stop("`time_labels` must have one entry per time point.", call. = FALSE)

  storage.mode(X) <- "double"
  colnames(X) <- .longitudinal_colnames(item_names, time_labels)

  # Waves with no observed item at all. These contribute a flat likelihood at
  # that occasion and are carried through by FIML rather than deleted.
  wave_missing <- matrix(FALSE, nrow = nrow(X), ncol = Tn)
  for (t in seq_len(Tn)) {
    block <- X[, ((t - 1) * J + 1):(t * J), drop = FALSE]
    wave_missing[, t] <- rowSums(!is.na(block)) == 0L
  }

  list(
    X            = X,
    n_items      = J,
    n_times      = Tn,
    item_names   = item_names,
    time_labels  = time_labels,
    wave_missing = wave_missing,
    any_missing  = anyNA(X)
  )
}

# Column positions of the t-th time block in the canonical layout.
.time_block_cols <- function(t, n_items) {
  ((t - 1L) * n_items + 1L):(t * n_items)
}
