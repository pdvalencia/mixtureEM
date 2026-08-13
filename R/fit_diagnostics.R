# ==============================================================================
# Absolute fit and local-fit diagnostics
# ==============================================================================
#
# Three statistics that applied latent class researchers lean on, and that a
# categorical mixture model can produce almost for free:
#
#   absolute_fit()        global fit of the model to the response-pattern table
#   bivariate_residuals() local fit, pair of items by pair of items
#   classification_table() how much a modal assignment costs in accuracy
#
# The first two apply only to fully categorical measurement models, because
# both are statements about a contingency table. The third applies to any
# model, since it only reads the posterior.
#
# None of the three needs the full response-pattern table enumerated. That is
# worth stating explicitly, because the obvious implementation of the Pearson
# X^2 does enumerate it and becomes unusable at a dozen items:
#
#   X^2 = sum_cells (n_c - e_c)^2 / e_c
#       = sum_observed (n_c - e_c)^2 / e_c  +  sum_unobserved e_c
#       = sum_observed (n_c - e_c)^2 / e_c  +  (N - sum_observed e_c),
#
# since an unobserved cell contributes (0 - e_c)^2 / e_c = e_c and the model
# probabilities sum to one over the whole table. G^2 and the Cressie-Read
# statistic weight each cell by its *observed* count, so unobserved cells drop
# out of those two outright. And the bivariate residuals need only the
# two-way margin, which conditional independence gives in closed form:
#
#   P(y_a = r, y_b = s) = sum_k gamma_k p_a(r | k) p_b(s | k).
#
# So every quantity here is computed from the observed patterns and the
# parameters, and none of them scales with the size of the table.

# ------------------------------------------------------------------------------
# Reading category probabilities out of an arbitrary measurement model
# ------------------------------------------------------------------------------

# One entry per column of the model's data matrix, in data-column order; each
# entry holds the K x R matrix of response probabilities for that indicator and
# the category codes those columns correspond to. Returns NULL if any indicator
# is not categorical, which is the signal to the callers that a contingency
# table does not exist for this model.
#
# The category codes matter: Bernoulli indicators are stored as one probability
# per item and are coded 0/1 in the data, whereas multinoulli indicators are
# stored one-hot and coded 1..max_val. Carrying the codes alongside the
# probabilities is what lets the observed table be cross-tabulated against the
# expected one without the caller knowing which family it is looking at.
.categorical_item_probs <- function(mm) {
  if (inherits(mm, c("nested", "blocks"))) {
    subs <- lapply(mm$models, .categorical_item_probs)
    if (any(vapply(subs, is.null, logical(1)))) return(NULL)
    return(unlist(subs, recursive = FALSE, use.names = FALSE))
  }

  pis <- mm$parameters$pis
  if (is.null(pis)) return(NULL)      # gaussian, poisson, structured normal, ...

  M <- mm$max_val
  if (is.null(M)) {
    # Bernoulli: the stored parameter is P(y = 1), so the second column of the
    # returned matrix is the stored one and the first is its complement.
    return(lapply(seq_len(ncol(pis)), function(j)
      list(probs = cbind(1 - pis[, j], pis[, j]), categories = c(0, 1))))
  }

  lapply(seq_len(ncol(pis) %/% M), function(j)
    list(probs = pis[, ((j - 1L) * M + 1L):(j * M), drop = FALSE],
         categories = seq_len(M)))
}

# The categorical items of a fitted model, aligned to its data columns and
# named after them, or NULL. `n_cols` guards against a measurement model whose
# parameter layout has drifted from the data it was fitted to.
.fit_item_probs <- function(mm, n_cols, col_names = NULL) {
  items <- .categorical_item_probs(mm)
  if (is.null(items) || length(items) != n_cols) return(NULL)
  nms <- col_names %||% paste0("Item", seq_len(n_cols))
  for (j in seq_len(n_cols)) items[[j]]$name <- nms[j]
  items
}

# Marginal class proportions. A covariate structural model makes class
# membership case-specific, at which point there is no single gamma to build a
# model-implied table from; the callers refuse those models rather than
# quietly averaging over cases.
.marginal_class_weights <- function(object) {
  if (has_covariate(object$sm)) return(NULL)
  object$weights
}

# ------------------------------------------------------------------------------
# Absolute fit
# ------------------------------------------------------------------------------

#' Absolute Fit of a Categorical Mixture Model
#'
#' @description
#' Compares the observed response-pattern frequencies with those the model
#' implies, over the contingency table formed by crossing every categorical
#' indicator. Three members of the power-divergence family are reported:
#' the likelihood-ratio statistic \eqn{G^2} (also written \eqn{L^2}), the
#' Pearson \eqn{X^2}, and the Cressie-Read statistic (\eqn{\lambda = 2/3}),
#' each on \eqn{df = W - P - 1} degrees of freedom, where \eqn{W} is the number
#' of cells in the table and \eqn{P} the number of free parameters.
#'
#' The statistics are defined only for fully categorical indicators observed
#' without missingness. Even then they should be read with care: the table has
#' \eqn{W} cells and is usually extremely sparse, so the chi-square reference
#' distribution is unreliable and the value is best used to compare models
#' rather than to test one in isolation. [`blrt()`] tests a model against one
#' with fewer classes without relying on that reference distribution.
#'
#' @param object A model fitted by [`fit_mixture()`], [`fit_rmlca()`] or
#'   [`fit_lta()`].
#' @return An object of class `absolute_fit` with elements `g2`, `x2`,
#'   `cressie_read`, `df`, the corresponding `p_value`s, `n_cells` and
#'   `n_patterns`; or `NULL` (with a message) when the statistics do not apply.
#' @seealso [`bivariate_residuals()`] for the local counterpart,
#'   [`classification_table()`], [`blrt()`].
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' absolute_fit(fit)
#'
#' @references
#' Langeheine, R., Pannekoek, J., & van de Pol, F. (1996). Bootstrapping
#' goodness-of-fit measures in categorical data analysis. \emph{Sociological
#' Methods & Research}, \emph{24}(4), 492-516.
#' \doi{10.1177/0049124196024004004}
#' @export
absolute_fit <- function(object) {
  info <- .nested_fit_info(object)
  X    <- info$data
  if (is.null(X)) {
    message("The fitted object does not retain its data.")
    return(NULL)
  }
  if (anyNA(X)) {
    message("Absolute fit is not defined with missing data; the contingency ",
            "table is incomplete. Compare models with BIC or a ",
            "likelihood-ratio test instead.")
    return(NULL)
  }
  if (isTRUE(info$conditional)) {
    message("Absolute fit is not defined once covariates enter the model: a ",
            "case's implied probability depends on its covariates, so there ",
            "is no single model-implied table to compare with. Refit the ",
            "measurement model on its own, or use `n_steps = 3`.")
    return(NULL)
  }

  levels_per_col <- .longitudinal_col_levels(info$mm, ncol(X))
  if (is.null(levels_per_col)) {
    message("Absolute fit requires categorical indicators throughout.")
    return(NULL)
  }

  w   <- info$weights
  key <- apply(X, 1, paste, collapse = "\r")
  obs <- tapply(w, key, sum)
  # One representative case per pattern carries that pattern's log-likelihood.
  rep_idx <- !duplicated(key)
  p_hat   <- exp(info$ll_case[rep_idx])
  names(p_hat) <- key[rep_idx]
  p_hat   <- p_hat[names(obs)]

  N   <- sum(w)
  exp_counts <- N * p_hat

  g2 <- 2 * sum(obs * log(obs / exp_counts))
  # The cells with no observed case contribute (0 - e)^2 / e = e apiece, and
  # they are exactly the mass the observed cells leave over.
  x2 <- sum((obs - exp_counts)^2 / exp_counts) + max(N - sum(exp_counts), 0)
  # Cressie-Read at lambda = 2/3, the read-across between the other two. Every
  # term carries the observed count as a factor, so unobserved cells vanish.
  lambda <- 2 / 3
  cr <- 2 / (lambda * (lambda + 1)) * sum(obs * ((obs / exp_counts)^lambda - 1))

  W  <- prod(levels_per_col)
  df <- W - info$n_params - 1
  p_of <- function(stat)
    if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_

  structure(
    list(g2 = g2, x2 = x2, cressie_read = cr, df = df,
         p_value = p_of(g2), p_value_x2 = p_of(x2),
         p_value_cressie_read = p_of(cr),
         n_cells = W, n_patterns = length(obs), n_params = info$n_params,
         label = info$label),
    class = "absolute_fit")
}

#' @export
print.absolute_fit <- function(x, ...) {
  cat("=========================================================\n")
  cat("                  ABSOLUTE FIT                           \n")
  cat("=========================================================\n")
  cat(sprintf("Table: %d cells, %d observed response patterns\n",
              x$n_cells, x$n_patterns))
  cat(sprintf("Free parameters: %d   df: %d\n\n", x$n_params, x$df))
  cat(sprintf("%-16s %12s %10s\n", "Statistic", "Value", "p-value"))
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  row <- function(lab, v, p)
    cat(sprintf("%-16s %12.4f %10s\n", lab, v,
                if (is.na(p)) "--" else sprintf("%.4f", p)))
  row("L-squared", x$g2, x$p_value)
  row("X-squared", x$x2, x$p_value_x2)
  row("Cressie-Read", x$cressie_read, x$p_value_cressie_read)
  cat("=========================================================\n")
  if (x$n_patterns < x$n_cells / 2)
    cat("Note: the table is sparse, so the chi-square reference\n",
        "distribution is unreliable. Prefer these statistics for\n",
        "comparing models over testing one in isolation.\n")
  invisible(x)
}

# ------------------------------------------------------------------------------
# Bivariate residuals
# ------------------------------------------------------------------------------

#' Bivariate Residuals
#'
#' @description
#' A local measure of fit: for each pair of categorical indicators, the Pearson
#' chi-square of the observed two-way table against the table the model implies,
#' divided by its degrees of freedom \eqn{(R_a - 1)(R_b - 1)}. If the model
#' were true, a bivariate residual should not be substantially larger than 1.
#'
#' Where [`absolute_fit()`] says whether the model fits, this says *where* it
#' fails: a large value flags the specific pair of items whose association the
#' latent classes do not reproduce, which is the conditional-independence
#' assumption showing its seams. This is the classic local-dependence
#' diagnostic (Oberski et al., 2013), and what the usual
#' referee question about local dependence asks for.
#'
#' Unlike the absolute-fit statistics, bivariate residuals do not require the
#' full response-pattern table and so remain usable with many indicators: the
#' model-implied two-way margin follows in closed form from conditional
#' independence, \eqn{P(y_a = r, y_b = s) = \sum_k \gamma_k\, p_a(r|k)\,
#' p_b(s|k)}.
#'
#' With missing data each pair is computed on the cases observing both items,
#' and the expected counts are scaled to that pair's total. This is a
#' pairwise-complete statistic rather than a full-information one, so read it
#' as descriptive when missingness is heavy.
#'
#' @param object A model fitted by [`fit_mixture()`] or [`fit_rmlca()`] with
#'   categorical indicators.
#' @return An object of class `bivariate_residuals`: a lower-triangular
#'   indicator-by-indicator matrix, `NA` on and above the diagonal; or `NULL`
#'   (with a message) when the statistic does not apply.
#' @seealso [`absolute_fit()`], [`classification_table()`].
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' bivariate_residuals(fit)
#'
#' @references
#' Oberski, D. L., van Kollenburg, G. H., & Vermunt, J. K. (2013). A Monte
#' Carlo evaluation of three methods to detect local dependence in binary
#' data latent class models. \emph{Advances in Data Analysis and
#' Classification}, \emph{7}(3), 267-279. \doi{10.1007/s11634-013-0146-2}
#' @export
bivariate_residuals <- function(object) {
  if (inherits(object, "lta_model")) {
    message("Bivariate residuals are not available for latent transition ",
            "models: the two-way margin of a pair of items at different ",
            "occasions runs through the transition matrices rather than a ",
            "single class distribution. Use `fit_rmlca()` for a ",
            "repeated-measures model whose residuals are defined.")
    return(NULL)
  }
  if (!inherits(object, "mixture_model"))
    stop("`object` must be a fitted mixture model.", call. = FALSE)

  X <- object$data
  if (is.null(X)) {
    message("The fitted object does not retain its data.")
    return(NULL)
  }
  gamma <- .marginal_class_weights(object)
  if (is.null(gamma)) {
    message("Bivariate residuals are not defined once covariates enter the ",
            "model: class membership is case-specific, so there is no single ",
            "model-implied two-way table. Refit the measurement model on its ",
            "own, or use `n_steps = 3`.")
    return(NULL)
  }

  items <- .fit_item_probs(object$mm, ncol(X), colnames(X))
  if (is.null(items)) {
    message("Bivariate residuals require categorical indicators throughout.")
    return(NULL)
  }

  J   <- length(items)
  w   <- object$sample_weights
  nms <- vapply(items, `[[`, character(1), "name")
  bvr <- matrix(NA_real_, J, J, dimnames = list(nms, nms))

  for (b in seq_len(J)) for (a in seq_len(b - 1L)) {
    ia <- items[[a]]; ib <- items[[b]]
    keep <- !is.na(X[, a]) & !is.na(X[, b])
    if (!any(keep)) next

    # Observed weighted two-way table, laid out over the model's own category
    # codes so that a category no case happens to take still gets a cell.
    fa <- factor(X[keep, a], levels = ia$categories)
    fb <- factor(X[keep, b], levels = ib$categories)
    obs <- tapply(w[keep], list(fa, fb), sum)
    obs[is.na(obs)] <- 0

    # Model-implied margin: t(P_a) diag(gamma) P_b, which is the closed form of
    # sum_k gamma_k p_a(r|k) p_b(s|k).
    prob <- crossprod(ia$probs * gamma, ib$probs)
    expected <- sum(w[keep]) * prob

    # A cell the model calls impossible contributes nothing if no case took it,
    # and an infinite chi-square if one did. Both are the right answer; what is
    # not is the NaN that 0/0 would otherwise put in place of the first, taking
    # the whole residual with it.
    zero <- expected == 0
    chisq <- sum((obs[!zero] - expected[!zero])^2 / expected[!zero]) +
      if (any(zero & obs > 0)) Inf else 0
    bvr[b, a] <- chisq / ((length(ia$categories) - 1L) *
                            (length(ib$categories) - 1L))
  }

  class(bvr) <- c("bivariate_residuals", "matrix", "array")
  bvr
}

#' @export
print.bivariate_residuals <- function(x, digits = 4, ...) {
  m <- unclass(x)
  J <- nrow(m)
  cat("=========================================================\n")
  cat("               BIVARIATE RESIDUALS                       \n")
  cat("=========================================================\n")
  cat("Pearson chi-square per item pair, divided by its df.\n")
  cat("Values well above 1 flag a pair whose association the\n")
  cat("classes do not reproduce (local dependence).\n\n")

  if (J < 2L) {
    cat("A single indicator: no pairs.\n")
    cat("=========================================================\n")
    return(invisible(x))
  }

  # Lower triangle only: the matrix is symmetric by
  # construction and the diagonal is not a residual.
  lab_w <- max(nchar(rownames(m)))
  col_w <- max(9L, max(nchar(colnames(m)[-J])) + 1L)
  cat(sprintf("%-*s", lab_w, ""))
  for (j in seq_len(J - 1L)) cat(sprintf("%*s", col_w, colnames(m)[j]))
  cat("\n")
  for (i in 2:J) {
    cat(sprintf("%-*s", lab_w, rownames(m)[i]))
    for (j in seq_len(i - 1L))
      cat(sprintf("%*s", col_w, formatC(m[i, j], format = "f", digits = digits)))
    cat("\n")
  }

  worst <- which(m == max(m, na.rm = TRUE), arr.ind = TRUE)[1, ]
  cat(sprintf("\nLargest: %s x %s = %.4f\n",
              rownames(m)[worst[1]], colnames(m)[worst[2]],
              m[worst[1], worst[2]]))
  cat("=========================================================\n")
  invisible(x)
}

#' Plot Bivariate Residuals as a Heat Table
#'
#' @description
#' The same lower triangle the print method shows, shaded so that the pairs
#' needing attention are visible at a glance rather than read off a table.
#'
#' The colour scale is anchored at 1, the value a residual should not much
#' exceed if the model is true, rather than at the largest residual present: a
#' scale stretched to fit the data would make a well-fitting model look exactly
#' like a badly fitting one. Cells at or below 1 are pale, cells above it
#' redden, and everything past `max_shade` shares the deepest colour so that a
#' single extreme pair cannot flatten the rest.
#'
#' @param x An object returned by [`bivariate_residuals()`].
#' @param max_shade Value at which the colour scale saturates.
#' @param main Plot title.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' plot(bivariate_residuals(fit))
#' @export
plot.bivariate_residuals <- function(x, max_shade = 4, main = NULL, ...) {
  m <- unclass(x)
  J <- nrow(m)
  if (J < 2L) {
    message("A single indicator: no pairs to plot.")
    return(invisible(x))
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(2, max(4, 0.55 * max(nchar(rownames(m))) + 1), 5, 2))

  # Drop the empty first column and last row so the triangle fills the panel.
  tri <- m[-1L, -J, drop = FALSE]
  nr  <- nrow(tri); nc <- ncol(tri)

  shades <- grDevices::colorRampPalette(c("#F7F7F7", "#FDDBC7", "#D6604D",
                                          "#67001F"))(64)
  # image() draws the first row at the bottom, so reverse to make the table
  # read in the same order as it prints; and clamp rather than rely on zlim,
  # which drops out-of-range cells entirely and would leave the worst pair in
  # the figure blank.
  shown <- pmin(tri, max_shade)
  graphics::image(seq_len(nc), seq_len(nr), t(shown[nr:1, , drop = FALSE]),
                  col = shades, zlim = c(0, max_shade), axes = FALSE,
                  xlab = "", ylab = "",
                  main = main %||% "Bivariate residuals")
  graphics::axis(3, at = seq_len(nc), labels = colnames(tri), tick = FALSE,
                 line = -0.5, cex.axis = 0.9)
  graphics::axis(2, at = seq_len(nr), labels = rev(rownames(tri)),
                 tick = FALSE, las = 1, line = -0.5, cex.axis = 0.9)

  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    v <- tri[i, j]
    if (is.na(v)) next
    graphics::text(j, nr - i + 1, formatC(v, format = "f", digits = 2),
                   col = if (v > max_shade / 2) "white" else "grey20",
                   cex = 0.9)
  }
  graphics::box()
  graphics::mtext(sprintf("shaded from 0 to %g; values above 1 flag local dependence",
                          max_shade),
                  side = 1, line = 0.5, cex = 0.8, col = "grey35")
  invisible(x)
}

# ------------------------------------------------------------------------------
# Classification table
# ------------------------------------------------------------------------------

#' Classification Table and Classification Error
#'
#' @description
#' Cross-classifies the probabilistic class memberships against the modal
#' assignment, which is what quantifies the cost of treating a fitted class as
#' though it were an observed group. Entry \eqn{(k, m)} is
#' \eqn{\sum_{i:\,\mathrm{modal}(i) = m} w_i P(k \mid y_i)}, so the rows sum to
#' the model-expected class sizes and the columns to the modal counts. The two
#' sets of totals disagree by exactly the amount modal assignment distorts the
#' class proportions.
#'
#' The classification error is \eqn{1 - \mathrm{trace}/N}: the proportion of
#' cases the modal rule is expected to place in the wrong class. It is the
#' quantity that motivates the bias-adjusted 3-step estimators, since it is
#' the error those corrections exist to undo — see [`fit_mixture()`]'s
#' `n_steps` and `correction` arguments.
#'
#' @param object A model fitted by [`fit_mixture()`], or by [`fit_lta()`], in
#'   which case one table per occasion is returned.
#' @return An object of class `classification_table`: the K x K matrix, with
#'   the classification `error`, `n` and expected/modal class sizes attached as
#'   attributes. For a latent transition model, a list of such tables.
#' @seealso [`classification_diagnostics()`], which prints this alongside the
#'   average posterior probabilities; [`class_assignments()`] for the per-case
#'   assignment this table aggregates.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' classification_table(fit)
#' @export
classification_table <- function(object) {
  if (inherits(object, "lta_model")) {
    tabs <- lapply(seq_len(object$longitudinal$n_times), function(t)
      .classification_table(object$gamma[[t]], object$weights_vec,
                            object$n_statuses))
    names(tabs) <- object$longitudinal$time_labels
    return(tabs)
  }
  if (!inherits(object, "mixture_model"))
    stop("`object` must be a fitted mixture model.", call. = FALSE)
  .classification_table(exp(object$log_resp), object$sample_weights,
                        object$n_components)
}

# resp: n x K posterior probabilities; w: case weights.
.classification_table <- function(resp, w, K) {
  modal <- max.col(resp, ties.method = "first")
  tab   <- matrix(0, K, K)
  for (m in seq_len(K)) {
    idx <- modal == m
    if (any(idx)) tab[, m] <- colSums(resp[idx, , drop = FALSE] * w[idx])
  }
  dimnames(tab) <- list(paste("Class", seq_len(K)),
                        paste("Modal", seq_len(K)))
  N <- sum(w)
  structure(tab, class = c("classification_table", "matrix", "array"),
            error = 1 - sum(diag(tab)) / N, n = N,
            expected_size = rowSums(tab), modal_size = colSums(tab))
}

#' @export
print.classification_table <- function(x, ...) {
  tab <- unclass(x)
  attr(tab, "error") <- attr(tab, "n") <- NULL
  attr(tab, "expected_size") <- attr(tab, "modal_size") <- NULL
  N <- attr(x, "n")

  cat("=========================================================\n")
  cat("               CLASSIFICATION TABLE                      \n")
  cat("=========================================================\n")
  cat("Rows: model-expected membership | Columns: modal assignment\n\n")
  full <- cbind(tab, Total = attr(x, "expected_size"))
  full <- rbind(full, Total = c(attr(x, "modal_size"), N))
  print(round(full, 4))
  cat(sprintf("\nClassification error: %.4f (%.2f%% of %.4g cases)\n",
              attr(x, "error"), 100 * attr(x, "error"), N))
  cat("=========================================================\n")
  invisible(x)
}
