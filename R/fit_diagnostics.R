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
#' The statistics require fully categorical indicators. Even then they should
#' be read with care: the table has \eqn{W} cells and is usually extremely
#' sparse, so the chi-square reference distribution is unreliable and the
#' value is best used to compare models rather than to test one in isolation.
#' [`blrt()`] tests a model against one with fewer classes without relying on
#' that reference distribution.
#'
#' @section Missing data:
#' With one or more missing values (categorical, plain `fit_mixture()` models
#' only), the statistics are computed under the missing-at-random (MAR)
#' assumption instead: the model is compared not to the raw response table,
#' which no longer exists once cases have different items observed, but to a
#' saturated model fit to the same partition of the data by which items each
#' case observed. `df` is smaller than in the complete-data case
#' (\eqn{df = W - 1 - P}) because the saturated baseline already accounts for
#' the missingness pattern. A short block giving the model's fit jointly with
#' the stronger missing-completely-at-random (MCAR) assumption is printed
#' underneath; use [`mcar_test()`] to test that assumption on its own. This
#' can be slow, or refused outright, once the number of indicator categories
#' crossed together grows large -- the same \eqn{W} that already makes the
#' complete-data table sparse.
#'
#' @param object A model fitted by [`fit_mixture()`], [`fit_rmlca()`] or
#'   [`fit_lta()`].
#' @return An object of class `absolute_fit` with elements `g2`, `x2`,
#'   `cressie_read`, `df`, the corresponding `p_value`s, `dissimilarity`,
#'   `n_cells` and `n_patterns`; or `NULL` (with a message) when the
#'   statistics do not apply. With missing data, also `g2_mcar`, `x2_mcar`,
#'   `cressie_read_mcar`, `df_mcar` and `p_value_mcar` for the block computed
#'   under MCAR, `ll_sat` for the saturated model's log-likelihood, and
#'   `mar = TRUE`.
#' @seealso [`bivariate_residuals()`] for the local counterpart,
#'   [`classification_table()`], [`blrt()`], [`mcar_test()`] for testing the
#'   missingness mechanism on its own.
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
  if (isTRUE(info$conditional)) {
    message("Absolute fit is not defined once covariates enter the model: a ",
            "case's implied probability depends on its covariates, so there ",
            "is no single model-implied table to compare with. Refit the ",
            "measurement model on its own, or use `n_steps = 3`.")
    return(NULL)
  }

  levels_per_col <- .longitudinal_col_levels(info$mm, ncol(X))
  if (is.null(levels_per_col)) {
    message("Absolute fit requires categorical indicators throughout: there ",
            "is no response-pattern contingency table to compare against ",
            "once an indicator is continuous. See `bivariate_residuals()` ",
            "for a local diagnostic that is defined for continuous ",
            "indicators, or `blrt()` for a global test that does not need ",
            "one.")
    return(NULL)
  }

  if (anyNA(X)) return(.absolute_fit_mar(object, info, X, levels_per_col))

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
  # Dissimilarity index: half the total absolute deviation between observed
  # and expected counts, expressed as a share of N. The unobserved cells
  # contribute their expected mass apiece, the same shortcut X^2 uses above.
  di <- (sum(abs(obs - exp_counts)) + max(N - sum(exp_counts), 0)) / (2 * N)

  W  <- prod(levels_per_col)
  df <- W - info$n_params - 1
  p_of <- function(stat)
    if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_

  structure(
    list(g2 = g2, x2 = x2, cressie_read = cr, df = df,
         p_value = p_of(g2), p_value_x2 = p_of(x2),
         p_value_cressie_read = p_of(cr), dissimilarity = di,
         n_cells = W, n_patterns = length(obs), n_params = info$n_params,
         label = info$label),
    class = "absolute_fit")
}

# ------------------------------------------------------------------------------
# Absolute fit with missing data (MAR)
# ------------------------------------------------------------------------------
#
# Everything here is scored on the missingness-pattern partition: cases are
# grouped by which items they observed, and within each group the cells are
# every combination those observed items can take, not just the ones anyone
# actually gave. Two probability models are scored against that same
# partition -- the fitted mixture, and a saturated model with one free
# probability per cell of the full W-cell table, fit by its own small EM
# (`.saturated_mar()`) so it best explains the missingness-pattern-partitioned
# data without assuming a mixture structure. The mixture-vs-observed row is
# the MCAR statistics (the model tested jointly with MCAR); the
# saturated-vs-observed row is a test of MCAR alone; and the MAR statistics
# this function heads with are the difference of the two, on
# df = W - 1 - n_params, which is the actual test of the model once the
# missingness mechanism is no longer part of what is being judged.

# Fits the saturated model under MAR: one probability per cell of the full,
# completely-crossed W-cell table, estimated by EM on the (missingness
# pattern, observed responses) partition of the data. Returns NULL, the same
# signal `absolute_fit()` already gives up on, when the full table is too
# large to enumerate.
.saturated_mar <- function(X, levels, w) {
  W <- prod(vapply(levels, length, integer(1)))
  n_patterns <- nrow(unique(is.na(X)))
  if (W > 1e5 || W * n_patterns > 2e7) return(NULL)

  cells <- as.matrix(expand.grid(levels))      # W x J, any consistent order
  key   <- apply(X, 1, function(r)
             paste(ifelse(is.na(r), ".", r), collapse = "|"))
  np    <- tapply(w, key, sum); pats <- names(np); np <- as.numeric(np)
  comp  <- matrix(TRUE, length(pats), nrow(cells))
  for (p in seq_along(pats)) {
    r <- X[match(pats[p], key), ]
    for (j in seq_len(ncol(X)))
      if (!is.na(r[j])) comp[p, ] <- comp[p, ] & (cells[, j] == r[j])
  }
  N <- sum(np); pi <- rep(1 / nrow(cells), nrow(cells)); old <- -Inf
  for (it in 1:5000) {
    num <- sweep(comp, 2, pi, "*"); den <- rowSums(num)
    ll  <- sum(np * log(den))
    pi  <- colSums(sweep(num / den, 1, np, "*")) / N
    if (abs(ll - old) < 1e-12) break
    old <- ll
  }
  list(pi = pi, ll = ll, cells = cells)
}

# All combinations a subset of columns can take, in each item's own category
# codes, so the result cross-tabulates against the raw data without the
# caller knowing which family (Bernoulli, multinoulli, ...) it is looking at.
.enumerate_cells <- function(items, cols)
  as.matrix(expand.grid(lapply(items[cols], `[[`, "categories")))

# P(cell) under the fitted mixture, for a subset of columns: conditional
# independence lets the unobserved columns simply not appear in the product,
# which is what makes this well-defined without integrating anything out.
.model_cell_prob <- function(items, gamma, cols, cells) {
  prob_k <- matrix(1, nrow(cells), length(gamma))
  for (jj in seq_along(cols)) {
    j   <- cols[jj]
    rho <- items[[j]]$probs
    ix  <- match(cells[, jj], items[[j]]$categories)
    prob_k <- prob_k * t(rho[, ix, drop = FALSE])
  }
  as.vector(prob_k %*% gamma)
}

# P(cell) under the saturated model, for a subset of columns: the full W-cell
# table is grouped by its values on that subset and summed, which marginalises
# out the columns not in the subset.
.sat_cell_prob <- function(full_cells, pi, cols, cells) {
  key_full  <- apply(full_cells[, cols, drop = FALSE], 1, paste, collapse = "|")
  key_cells <- apply(cells, 1, paste, collapse = "|")
  unname(tapply(pi, key_full, sum)[key_cells])
}

# The four statistics, accumulated cell by cell over every missingness
# pattern's own sub-table, for both the model and the saturated baseline.
.mar_partition_stats <- function(X, items, gamma, sat, w) {
  miss   <- is.na(X)
  patkey <- apply(miss, 1, function(r) paste(as.integer(r), collapse = ""))
  upat   <- unique(patkey)
  lambda <- 2 / 3

  g2_m <- x2_m <- cr_m <- di_m <- 0
  g2_s <- x2_s <- cr_s <- 0
  s_terms <- 0L

  for (pk in upat) {
    idx <- which(patkey == pk)
    obs_cols <- which(!miss[idx[1L], ])
    if (length(obs_cols) == 0L) next   # nothing observed, nothing to score
    n_m <- sum(w[idx])

    cells   <- .enumerate_cells(items, obs_cols)
    subkey  <- apply(X[idx, obs_cols, drop = FALSE], 1, paste, collapse = "|")
    cellkey <- apply(cells, 1, paste, collapse = "|")
    o <- tapply(w[idx], factor(subkey, levels = cellkey), sum)
    o <- ifelse(is.na(o), 0, as.numeric(o))

    e_model <- n_m * .model_cell_prob(items, gamma, obs_cols, cells)
    e_sat   <- n_m * .sat_cell_prob(sat$cells, sat$pi, obs_cols, cells)

    # As in the complete-data statistics: a cell no one took (o == 0)
    # contributes to X^2 and DI but not to L2 or Cressie-Read. A cell the
    # saturated model calls impossible (e_sat == 0) is skipped from X^2 rather
    # than counted as an infinite or undefined term -- the saturated model has
    # one free probability per full-table cell, so a cell can end up on the
    # boundary the same way a fitted class probability can.
    nz  <- o > 0
    nzm <- e_model > 0
    nzs <- e_sat > 0
    g2_m <- g2_m + 2 * sum(o[nz] * log(o[nz] / e_model[nz]))
    g2_s <- g2_s + 2 * sum(o[nz] * log(o[nz] / e_sat[nz]))
    x2_m <- x2_m + sum((o[nzm] - e_model[nzm])^2 / e_model[nzm])
    x2_s <- x2_s + sum((o[nzs] - e_sat[nzs])^2 / e_sat[nzs])
    cr_m <- cr_m + 2 / (lambda * (lambda + 1)) *
      sum(o[nz] * ((o[nz] / e_model[nz])^lambda - 1))
    cr_s <- cr_s + 2 / (lambda * (lambda + 1)) *
      sum(o[nz] * ((o[nz] / e_sat[nz])^lambda - 1))
    di_m <- di_m + sum(abs(o - e_model))

    s_terms <- s_terms + (nrow(cells) - 1L)
  }

  list(g2_model = g2_m, x2_model = x2_m, cr_model = cr_m,
       di_model = di_m / (2 * sum(w)),
       g2_sat = g2_s, x2_sat = x2_s, cr_sat = cr_s,
       s_terms = s_terms, n_patterns = length(upat))
}

# absolute_fit()'s branch for incomplete data. `object` is needed (rather than
# just `info`) because the marginal class weights are read off it directly,
# the same accessor bivariate_residuals() uses.
.absolute_fit_mar <- function(object, info, X, levels_per_col) {
  if (!inherits(object, "mixture_model")) {
    message("Absolute fit with missing data is available for a plain ",
            "fit_mixture() model only.")
    return(NULL)
  }
  gamma <- .marginal_class_weights(object)

  items <- .fit_item_probs(info$mm, ncol(X), colnames(X))
  w     <- info$weights
  p     <- info$n_params
  W     <- prod(levels_per_col)

  sat <- .saturated_mar(X, lapply(items, `[[`, "categories"), w)
  if (is.null(sat)) {
    message("Absolute fit with missing data needs the full crossing of ",
            "categories (", W, " cells here), which is too large to ",
            "enumerate. Compare models with BIC or a likelihood-ratio ",
            "test instead.")
    return(NULL)
  }
  st <- .mar_partition_stats(X, items, gamma, sat, w)

  df_mcar <- st$s_terms - p
  df_mar  <- W - 1L - p

  g2 <- st$g2_model - st$g2_sat
  x2 <- st$x2_model - st$x2_sat
  cr <- st$cr_model - st$cr_sat

  p_of <- function(stat, df)
    if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_

  structure(
    list(g2 = g2, x2 = x2, cressie_read = cr, df = df_mar,
         p_value = p_of(g2, df_mar), p_value_x2 = p_of(x2, df_mar),
         p_value_cressie_read = p_of(cr, df_mar),
         g2_mcar = st$g2_model, x2_mcar = st$x2_model,
         cressie_read_mcar = st$cr_model, df_mcar = df_mcar,
         p_value_mcar = p_of(st$g2_model, df_mcar),
         ll_sat = sat$ll, dissimilarity = st$di_model,
         n_cells = W, n_patterns = st$n_patterns, n_params = p,
         label = info$label, mar = TRUE),
    class = "absolute_fit")
}

#' @export
print.absolute_fit <- function(x, ...) {
  cat("=========================================================\n")
  cat("                  ABSOLUTE FIT                           \n")
  cat("=========================================================\n")
  if (isTRUE(x$mar))
    cat("Missing data: statistics computed under MAR, by comparing\n",
        "the model to a saturated baseline on the same partition\n",
        "of the data by missingness pattern. See `mcar_test()` to\n",
        "test the missingness mechanism itself.\n\n")
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
  if (!is.null(x$dissimilarity))
    cat(sprintf("%-16s %12.4f %10s\n", "Dissimilarity", x$dissimilarity, ""))
  cat("=========================================================\n")
  if (isTRUE(x$mar))
    cat(sprintf(
      "Model tested jointly with MCAR: L2 %.4f  X2 %.4f  CR %.4f  df %d\n",
      x$g2_mcar, x$x2_mcar, x$cressie_read_mcar, x$df_mcar))
  else if (x$n_patterns < x$n_cells / 2)
    cat("Note: the table is sparse, so the chi-square reference\n",
        "distribution is unreliable. Prefer these statistics for\n",
        "comparing models over testing one in isolation.\n")
  invisible(x)
}

#' Test Whether Data Are Missing Completely at Random
#'
#' @description
#' Tests the missingness mechanism itself, separately from whether the fitted
#' model fits. [`absolute_fit()`] already reports a model comparison under the
#' weaker missing-at-random (MAR) assumption; `mcar_test(fit)` asks the
#' stronger question a reviewer sometimes wants answered on its own -- whether
#' the pattern of missing values could plausibly be unrelated to the data,
#' rather than depending on it.
#'
#' The test compares the observed response frequencies, partitioned by which
#' items each case answered, against a saturated model fit to that same
#' partition. It does not depend on the mixture model fitting well: the
#' saturated model is as flexible as the data allow, so what remains is a
#' statement about the missingness pattern, not about the number of classes.
#'
#' @param object A model fitted by [`fit_mixture()`], with categorical
#'   indicators and at least one missing value.
#' @return A list with `stat`, `df` and `p_value`; or `NULL` (with a message)
#'   when the data are complete, since there is then nothing to test.
#' @section Reading the result:
#' A small `p_value` says the missingness is **not** missing completely at
#' random -- whether a value is missing depends on the data in some way. That
#' is a common and often unsurprising finding (people who skip a question
#' about drug use are not a random subset of respondents), and it does **not**
#' by itself mean the weaker, more common missing-at-random assumption fails
#' too; nothing here tests that. It also says nothing about whether the
#' mixture model itself fits -- that question belongs to
#' [`absolute_fit()`], which enters only through the size of the table this
#' test uses, not through its own log-likelihood.
#' @seealso [`absolute_fit()`], whose missing-data branch this function
#'   reuses.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), ncol = 6)
#' X[sample(length(X), 30)] <- NA
#' fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#' mcar_test(fit)
#' @export
mcar_test <- function(object) {
  info <- .nested_fit_info(object)
  X    <- info$data
  if (is.null(X) || !anyNA(X)) {
    message("The data are complete: there is no missingness to test.")
    return(NULL)
  }

  af <- absolute_fit(object)
  if (is.null(af) || !isTRUE(af$mar)) return(NULL)

  stat <- af$g2_mcar - af$g2
  df   <- af$df_mcar - af$df
  p    <- if (df > 0) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_

  structure(list(stat = stat, df = df, p_value = p), class = "mcar_test")
}

#' @export
print.mcar_test <- function(x, ...) {
  cat("=========================================================\n")
  cat("               MISSING COMPLETELY AT RANDOM              \n")
  cat("=========================================================\n")
  cat(sprintf("Chi-square: %.4f   df: %d   p-value: %s\n", x$stat, x$df,
              if (is.na(x$p_value)) "--" else sprintf("%.4f", x$p_value)))
  cat("=========================================================\n")
  if (!is.na(x$p_value) && x$p_value < 0.05)
    cat("Rejected: the missingness depends on the data in some way.\n",
        "This says nothing about whether MAR holds, and nothing about\n",
        "whether the mixture model itself fits.\n")
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
#' @section How much to trust it:
#' The chi-square reference for this statistic does not work, and the evidence
#' is blunt. Over the eight null conditions of Oberski et al. (2013, Table 1) --
#' loadings of .5 and .8 crossed with n of 200, 500, 1000 and 5000, 200 samples
#' each -- a bivariate residual referred to chi-square rejected at a nominal 5
#' percent level in **zero of 200 samples in every one of the eight**. Its
#' empirical mean ran between 0.25 and 0.36 against the 1 that a
#' \eqn{\chi^2_1} has, and its variance between 0.1 and 0.2 against 2. Three
#' consequences follow, and all three matter more than the usual hedging
#' suggests:
#'
#' \itemize{
#'   \item **A low bivariate residual is not evidence of good fit.** A statistic
#'     that never rejects when the model is true also has nothing to say when it
#'     is. This is the authors' own closing point.
#'   \item **The ranking is weaker than it looks.** In their Figure 1 the naive
#'     bivariate residual has uniformly the lowest power of the three methods
#'     compared. It is adequate only against a residual correlation of about
#'     -0.4; it needs n of 5000 or more for correlations of \eqn{\pm 0.2} and
#'     -0.2, and it almost never detects \eqn{\pm 0.05}. All three methods lose
#'     power as the loadings grow, so well-separated classes hide local
#'     dependence rather than expose it.
#'   \item **With missing data a large value is ambiguous.** Under MAR the
#'     observed side of any residual statistic carries selection bias
#'     (Asparouhov & Muthen, 2015), so a large residual may be reporting the
#'     missingness mechanism rather than local dependence.
#' }
#'
#' `n_reps` replaces the broken reference distribution with a parametric
#' bootstrap, which in the same simulation held between 0.020 and 0.085 against
#' a nominal 0.05. Use it before drawing any conclusion from the size of a
#' residual. Each bootstrap replicate is blanked out in the same cells as the
#' real data, so a p-value from missing data is exact when the data are
#' missing completely at random and an approximation otherwise -- the
#' replicate's missingness cannot reproduce a dependence on class or response
#' that the real gaps might carry.
#'
#' Note that the statistic bootstrapped is the one this function computes,
#' \eqn{\chi^2} divided by its degrees of freedom, and not a raw Pearson
#' \eqn{\chi^2}. For binary items the degrees of freedom are 1 and the two
#' coincide; for polytomous items they do not, so do not compare the number
#' printed here against a raw chi-square table. The bootstrap is applied to
#' whatever statistic is computed, so it is calibrated either way.
#'
#' With missing data each pair is computed on the cases observing both items,
#' and the expected counts are evaluated at those cases' own posterior class
#' membership rather than at the class proportions for the whole sample. This
#' matters whenever the two subsets differ -- an item that is missing more
#' often for one class than another otherwise makes an unrelated pair of items
#' look locally dependent, when what actually happened is that missingness
#' changed who is left in the comparison. This is a pairwise-complete
#' statistic rather than a full-information one, so still read it as
#' descriptive when missingness is heavy. For categorical indicators this is
#' the same statistic, with the same divisor, that another program reports as
#' a bivariate residual, and the values agree closely on the same fit.
#'
#' For a plain continuous (Gaussian) measurement model with no missing data, a
#' different statistic is returned instead: for each item pair and class, the
#' modification index (Sorbom, 1989) for freeing the within-class residual
#' covariance of that pair, with its variance adjusted for the model's other
#' parameters, following Oberski, van Kollenburg and Vermunt (2013). It is
#' computed separately in each class, because a residual dependence can run in
#' opposite directions in different classes and a pooled statistic would
#' average it away. Another program's "bivariate residual" for continuous
#' indicators is deliberately *not* adjusted for the model's other parameters,
#' so the two do not agree numerically; the modification index is preferred
#' here on the evidence of Oberski et al.'s simulation, in which a bivariate
#' residual referred to chi-square gave below-nominal size and inadequate
#' power, while the modification index reproduced its nominal distribution and
#' was the more powerful of the two adequate methods. Detection is reliable
#' only when the offending effects are few and the measurement model is
#' strong (Janssen, van Laar, de Rooij, Kuha & Bakk, 2019); unmodelled
#' within-class dependence is worth taking seriously mainly because it can
#' manufacture spurious classes (Bauer & Curran, 2004).
#'
#' @param object A model fitted by [`fit_mixture()`] or [`fit_rmlca()`] with
#'   categorical indicators, or with a plain continuous measurement model.
#' @param n_reps Number of parametric bootstrap replicates used to calibrate the
#'   residuals. The default, `0`, computes the residuals alone and is the
#'   cheaper, uncalibrated diagnostic. `100` is the recommended working value
#'   and gives a Monte Carlo standard error of about 0.022 at \eqn{p = 0.05},
#'   which is adequate for flagging a pair; `500`, the number Oberski et al.
#'   used, is the publication-grade setting. Ignored for a continuous
#'   measurement model.
#' @param n_init_boot Random starts per bootstrap replicate. Replicates are
#'   refit without the final refinement step, since a replicate needs a
#'   residual rather than polished estimates.
#' @param verbose Report bootstrap progress.
#' @return For categorical indicators, an object of class
#'   `bivariate_residuals`: a lower-triangular indicator-by-indicator matrix,
#'   `NA` on and above the diagonal. When `n_reps > 0` a matrix of bootstrap
#'   p-values, laid out the same way, is attached as the `"p"` attribute and
#'   printed beside each residual. The sum of all pairwise residuals is
#'   attached as the `"total"` attribute and printed as "Total BVR", a single
#'   headline number for how much local dependence the model as a whole is
#'   carrying. For a continuous measurement model, an
#'   object of class `bivariate_residuals_gaussian` holding the modification
#'   index and expected parameter change per pair per class, the model-implied
#'   residual covariance and correlation, and a count of pairs where the
#'   information matrix needed the outer-product fallback.
#'   `NULL` (with a message) when neither applies.
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
#'
#' Sorbom, D. (1989). Model modification. \emph{Psychometrika}, \emph{54},
#' 371-384.
#'
#' Janssen, J. H. M., van Laar, S., de Rooij, M., Kuha, J., & Bakk, Z. (2019).
#' The detection of local dependence in the presence of one continuous latent
#' variable: A comparison of different statistics. \emph{Structural Equation
#' Modeling}, \emph{26}(2), 280-290.
#'
#' Bauer, D. J., & Curran, P. J. (2004). The integration of continuous and
#' discrete latent variable models: Potential problems and promising
#' opportunities. \emph{Psychological Methods}, \emph{9}(1), 3-29.
#'
#' van Kollenburg, G. H., Mulder, J., & Vermunt, J. K. (2015). Assessing model
#' fit in latent class analysis when asymptotics do not hold.
#' \emph{Methodology}, \emph{11}(2), 65-79.
#'
#' Asparouhov, T., & Muthen, B. (2015). Residual associations in latent class
#' and latent transition analysis. \emph{Structural Equation Modeling},
#' \emph{22}(2), 169-177.
#' @export
bivariate_residuals <- function(object, n_reps = 0, n_init_boot = 10,
                                verbose = FALSE) {
  n_reps <- as.integer(n_reps)
  if (is.na(n_reps) || n_reps < 0L)
    stop("`n_reps` must be a non-negative number of bootstrap draws.",
         call. = FALSE)

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
    if (identical(class(object$mm)[1], "gaussian_diag")) {
      # The continuous branch returns a modification index, which already has a
      # usable reference distribution and needs no bootstrap.
      if (n_reps > 0L)
        message("`n_reps` is ignored for a continuous measurement model: the ",
                "modification index returned there reproduces its nominal ",
                "chi-square distribution, which is the reason the bootstrap ",
                "exists for the categorical statistic.")
      return(.bivariate_mi_gaussian(object))
    }
    message("Bivariate residuals require categorical indicators throughout, ",
            "or a plain continuous (Gaussian) measurement model with no ",
            "missing data, mixing with categorical items, or repeated-",
            "measures / growth structure.")
    return(NULL)
  }

  resp <- exp(object$log_resp)
  bvr <- .bvr_matrix(items, X, object$sample_weights, resp)

  if (n_reps > 0L)
    attr(bvr, "p") <- .bvr_bootstrap(object, bvr, n_reps, n_init_boot, verbose)

  attr(bvr, "total") <- sum(bvr, na.rm = TRUE)
  class(bvr) <- c("bivariate_residuals", "matrix", "array")
  bvr
}

# The residual matrix itself, given the pieces the caller has already resolved.
# Split out of bivariate_residuals() so that a bootstrap replicate is scored by
# exactly the same code as the observed data rather than by a copy of it.
#
# resp: n x K posterior class probabilities. The expected two-way margin is
# evaluated at the posterior composition of the cases actually retained for
# each pair, not at the marginal class weights gamma -- the two agree only
# when the retained subset is the whole sample, i.e. when neither item has
# missing data.
.bvr_matrix <- function(items, X, w, resp) {
  J   <- length(items)
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

    # Model-implied margin over the retained cases: t(P_a) diag(gk) P_b, where
    # gk = sum_i w_i P(class = k | y_i) is the posterior-weighted class total
    # among cases with both items observed.
    gk <- colSums(resp[keep, , drop = FALSE] * w[keep])
    expected <- crossprod(ia$probs * gk, ib$probs)

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

  bvr
}

# Parametric bootstrap p-values for the residual matrix. The reference
# distribution is generated the way blrt() generates its own: draw class
# memberships from the fitted class weights, generate responses from the fitted
# measurement model, refit, and score the replicate. The classes of a replicate
# need no alignment to the classes of the observed fit -- the residual is a
# function of the two-way margins, which are invariant to relabelling the
# classes, so align_classes() would only cost time.
.bvr_bootstrap <- function(object, obs, n_reps, n_init_boot, verbose) {
  K <- object$n_components
  N <- nrow(object$data)
  J <- nrow(obs)

  ge <- matrix(0L, J, J)   # replicates at least as extreme as the observed
  ok <- matrix(0L, J, J)   # replicates that produced a residual at all

  if (verbose)
    message(sprintf("Bivariate residuals: %d bootstrap draws...", n_reps))

  for (i in seq_len(n_reps)) {
    classes <- sample(seq_len(K), size = N, replace = TRUE,
                      prob = object$weights)
    X_gen   <- generate_synthetic_data(object$mm, classes, N)
    colnames(X_gen) <- colnames(object$data)
    # Match the observed missingness pattern, cell for cell, so the replicate
    # is scored on data shaped the same way as the observed fit. Exact under
    # MCAR; under MAR or MNAR this only approximates the true reference
    # distribution, since the replicate's missingness carries none of the
    # dependence on class or response that generated the real gaps.
    X_gen[is.na(object$data)] <- NA

    # refine = FALSE for the same reason blrt() uses it: a replicate needs a
    # residual, not polished estimates, and the refinement is the expensive part.
    rep_fit <- try(fit_mixture_internal(
      X = X_gen, n_components = K,
      measurement = object$measurement_descriptor,
      n_init = n_init_boot, refine = FALSE), silent = TRUE)
    if (inherits(rep_fit, "try-error")) next

    items_r <- .fit_item_probs(rep_fit$mm, ncol(X_gen), colnames(X_gen))
    gamma_r <- .marginal_class_weights(rep_fit)
    if (is.null(items_r) || is.null(gamma_r)) next

    boot <- .bvr_matrix(items_r, X_gen, rep_fit$sample_weights,
                        exp(rep_fit$log_resp))

    good <- !is.na(boot) & !is.na(obs)
    ok[good] <- ok[good] + 1L
    ge[good] <- ge[good] + (boot[good] >= obs[good])

    if (verbose && (i %% max(1L, n_reps %/% 10L) == 0L))
      message(sprintf("  %d / %d draws complete", i, n_reps))
  }

  p <- ifelse(ok > 0L, ge / ok, NA_real_)
  dimnames(p) <- dimnames(obs)
  p
}

#' @export
print.bivariate_residuals <- function(x, digits = 4, ...) {
  p <- attr(x, "p")
  m <- unclass(x)
  attr(m, "p") <- NULL
  J <- nrow(m)
  cat("=========================================================\n")
  cat("               BIVARIATE RESIDUALS                       \n")
  cat("=========================================================\n")
  cat("Pearson chi-square per item pair, divided by its df.\n")
  if (is.null(p)) {
    cat("Ranks which pairs strain the model. NOT a calibrated\n")
    cat("test: referred to chi-square this statistic almost\n")
    cat("never rejects, so a low value is not evidence of fit.\n")
    cat("Use n_reps for bootstrap p-values.\n\n")
  } else {
    cat("Bootstrap p-value in brackets (proportion of replicates\n")
    cat("at least as extreme). Small p flags local dependence.\n\n")
  }

  if (J < 2L) {
    cat("A single indicator: no pairs.\n")
    cat("=========================================================\n")
    return(invisible(x))
  }

  # Lower triangle only: the matrix is symmetric by
  # construction and the diagonal is not a residual.
  cell <- function(i, j) {
    v <- formatC(m[i, j], format = "f", digits = digits)
    if (is.null(p) || is.na(p[i, j])) v
    else sprintf("%s [%.2f]", v, p[i, j])
  }

  lab_w <- max(nchar(rownames(m)))
  col_w <- max(9L, max(nchar(colnames(m)[-J])) + 1L)
  if (!is.null(p)) col_w <- col_w + 7L   # room for the bracketed p-value
  cat(sprintf("%-*s", lab_w, ""))
  for (j in seq_len(J - 1L)) cat(sprintf("%*s", col_w, colnames(m)[j]))
  cat("\n")
  for (i in 2:J) {
    cat(sprintf("%-*s", lab_w, rownames(m)[i]))
    for (j in seq_len(i - 1L)) cat(sprintf("%*s", col_w, cell(i, j)))
    cat("\n")
  }

  worst <- which(m == max(m, na.rm = TRUE), arr.ind = TRUE)[1, ]
  cat(sprintf("\nLargest: %s x %s = %.4f\n",
              rownames(m)[worst[1]], colnames(m)[worst[2]],
              m[worst[1], worst[2]]))
  total <- attr(x, "total")
  if (!is.null(total)) cat(sprintf("Total BVR: %.4f\n", total))
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
