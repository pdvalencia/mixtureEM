# Merging fitted classes into clusters: merge_classes() (roadmap Part 31).
#
# Pure post-processing on the measurement model's posterior probabilities.
# Nothing is refitted: a merged solution is the same mixture with some of its
# components read together, so the likelihood, the information criteria and
# every parameter of the fit are those of the fit itself.

# x log x with the limit 0 log 0 = 0, so that an exact 0 or 1 posterior
# contributes nothing rather than NaN.
.xlogx <- function(p) {
  out <- p * log(p)
  out[p <= 0] <- 0
  out
}

# Baudry et al.'s (2010) greedy path on an n x G posterior matrix. At each step
# the two clusters whose union removes the most classification entropy are
# combined, their columns summed, and the search repeats on the K - 1 columns.
# The union of `a` and `b` takes column `a`'s place, so with a < b every
# cluster stays ordered by its smallest original class.
.merge_path <- function(resp, w = rep(1, nrow(resp))) {
  G      <- ncol(resp)
  groups <- as.list(seq_len(G))
  cur    <- resp
  ent    <- function(m) -sum(w * .xlogx(m))

  label  <- function(g) paste0("{", paste(g, collapse = ","), "}")
  post   <- vector("list", G)
  clus   <- vector("list", G)
  merged <- rep(NA_character_, G)
  moved  <- rep(NA_real_, G)
  entK   <- rep(NA_real_, G)

  record <- function(k) {
    labs <- vapply(groups, label, character(1))
    colnames(cur) <- labs
    post[[k]] <<- cur
    clus[[k]] <<- groups
    entK[k]   <<- ent(cur)
  }
  record(G)

  for (k in rev(seq_len(G))[-G]) {
    best <- c(NA, NA); gain <- -Inf
    for (a in seq_len(k)[-k]) for (b in (a + 1L):k) {
      s <- cur[, a] + cur[, b]
      g <- sum(w * (.xlogx(s) - .xlogx(cur[, a]) - .xlogx(cur[, b])))
      if (g > gain) { gain <- g; best <- c(a, b) }
    }
    a <- best[1]; b <- best[2]
    merged[k - 1L] <- paste(label(groups[[a]]), "+", label(groups[[b]]))
    moved[k - 1L]  <- sum(w * (cur[, a] + cur[, b]))
    cur[, a]    <- cur[, a] + cur[, b]
    cur         <- cur[, -b, drop = FALSE]
    groups[[a]] <- sort(c(groups[[a]], groups[[b]]))
    groups[[b]] <- NULL
    record(k - 1L)
  }

  path <- data.frame(
    K             = rev(seq_len(G)),
    merged        = rev(merged),
    merge_entropy = rev(entK),
    entropy_drop  = c(NA, -diff(rev(entK))),
    cases_merged  = rev(moved),
    stringsAsFactors = FALSE)
  list(path = path, clusters = rev(clus), posterior = rev(post))
}

#' Merge fitted classes into clusters
#'
#' Combines the classes of a fitted model, two at a time, into fewer clusters,
#' following the entropy criterion of Baudry et al. (2010), and returns the
#' whole path from the fitted number of classes down to one. Nothing is
#' refitted: each merged solution is the fitted mixture with some of its
#' classes read together, so it fits the data exactly as well as the fit
#' itself.
#'
#' @details
#' A Gaussian class is a symmetric bell. When a group of cases has some other
#' shape (skewed, flat-topped, with a long tail), the best-fitting mixture can
#' spend two or three classes describing that one group, and the information
#' criteria will prefer that larger model. Merging is one way to read such a
#' solution: keep the classes the criteria asked for, then combine the ones
#' that describe the same group. It needs only the posterior class
#' probabilities, so it applies to every measurement model the package fits.
#'
#' At each step the two clusters \eqn{k, k'} chosen are those maximising
#' \deqn{\sum_i w_i \{(t_{ik} + t_{ik'}) \log(t_{ik} + t_{ik'}) -
#'   t_{ik} \log t_{ik} - t_{ik'} \log t_{ik'}\},}
#' the drop in classification entropy from reading them as one; their
#' posterior columns are then added. The posteriors are those of the
#' measurement model: for a fit whose structural model was attached in a
#' later step they are recovered without it, as in
#' [`classification_diagnostics()`].
#'
#' **No number of clusters is chosen.** Which classes belong together is not
#' identified by the data (Hennig, 2010): it depends on what one takes a
#' cluster to be. Baudry et al. suggest looking for an elbow in the plot of the
#' entropy against K, or for the K where `entropy_drop` falls away, or taking
#' the K that ICL favours across separately fitted models (the `ICL` column of
#' [`compare_mixtures()`], printed with `columns = "all"`); all three are a
#' choice the analyst makes and reports.
#'
#' **This is not the entropy used against choosing the number of classes.**
#' The advice not to choose K by entropy concerns the *relative* entropy
#' reported in `fit$metrics$entropy`, an index of how well separated the
#' classes of competing models are. `merge_entropy` is the *absolute*
#' classification entropy of one fitted model read at different levels of
#' merging; the likelihood is the same at every row, and the column only
#' records how much ambiguity in the class assignments each merge resolves. At
#' the fitted K it equals \eqn{(\mathrm{ICL} - \mathrm{BIC}) / 2}.
#'
#' @param object A fitted \code{mixture_model} with at least two classes.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"class_merge"}: a list with
#'   \describe{
#'     \item{`path`}{A data frame with one row per number of clusters `K`,
#'       from the fitted number down to 1: `merged`, the two clusters combined
#'       to reach that row; `merge_entropy`, the absolute classification
#'       entropy at that K; `entropy_drop`, its fall from the row above; and
#'       `cases_merged`, the (weighted) number of cases carried by the two
#'       clusters combined.}
#'     \item{`clusters`}{A list, one element per row of `path`, giving the
#'       fitted classes that make up each cluster.}
#'     \item{`posterior`}{A list of the n x K posterior matrices of the merged
#'       solutions, columns named by their classes, e.g. `"{1,5}"`.}
#'   }
#'   The fitted model is not modified.
#'
#' @references
#' Baudry, J.-P., Raftery, A. E., Celeux, G., Lo, K., & Gottardo, R. (2010).
#' Combining mixture components for clustering. \emph{Journal of
#' Computational and Graphical Statistics}, \emph{19}(2), 332-353.
#' \doi{10.1198/jcgs.2010.08111}
#'
#' Hennig, C. (2010). Methods for merging Gaussian mixture components.
#' \emph{Advances in Data Analysis and Classification}, \emph{4}(1), 3-34.
#' \doi{10.1007/s11634-010-0058-3}
#'
#' @examples
#' set.seed(1)
#' X <- rbind(matrix(rnorm(200, 0), ncol = 2), matrix(rnorm(200, 4), ncol = 2))
#' fit <- fit_mixture(X, n_classes = 3, measurement = "continuous",
#'                    n_init = 2, random_state = 1)
#' m <- merge_classes(fit)
#' m
#' plot(m)
#'
#' @export
merge_classes <- function(object, ...) UseMethod("merge_classes")

#' @export
merge_classes.default <- function(object, ...) {
  if (is.null(object$log_resp))
    stop("`object` must be a fitted model from fit_mixture().", call. = FALSE)
  K <- object$n_components
  if (K < 2L)
    stop("A one-class model has nothing to merge.", call. = FALSE)
  resp <- .step1_resp(object)
  w    <- object$sample_weights %||% rep(1, nrow(resp))
  out  <- .merge_path(resp, w)
  class(out) <- "class_merge"
  out
}

#' @export
merge_classes.lta_model <- function(object, ...) {
  stop("merge_classes() does not merge latent statuses: combining two ",
       "statuses changes the transition probabilities, which the merged ",
       "posteriors alone do not give.", call. = FALSE)
}

#' @export
print.class_merge <- function(x, digits = 3, ...) {
  cat("Merge path (Baudry et al., 2010)\n")
  cat("The likelihood and every parameter are those of the fit at every K.\n\n")
  p <- x$path
  p$merged[is.na(p$merged)] <- "(fitted classes)"
  p$merge_entropy <- round(p$merge_entropy, digits)
  p$entropy_drop  <- round(p$entropy_drop, digits)
  p$cases_merged  <- round(p$cases_merged, 1)
  print(p, row.names = FALSE, na.print = "")
  cat("\nmerge_entropy is the absolute classification entropy, not the",
      "\nrelative entropy of print(fit). No K is chosen; see ?merge_classes.\n")
  invisible(x)
}

#' @export
plot.class_merge <- function(x, ...) {
  p <- x$path
  op <- par(mfrow = c(1, 2))
  on.exit(par(op))
  plot(p$K, p$merge_entropy, type = "b", pch = 19, xlab = "Number of clusters",
       ylab = "Classification entropy", main = "Entropy by K", xaxt = "n")
  axis(1, at = p$K)
  # Baudry et al.'s rescaled plot: entropy against the cumulative number of
  # cases merged, where a change of slope marks where merging stops being free.
  cum <- c(0, cumsum(p$cases_merged[-1]))
  plot(cum, p$merge_entropy, type = "b", pch = 19,
       xlab = "Cumulative cases merged", ylab = "Classification entropy",
       main = "Entropy by cases merged")
  text(cum, p$merge_entropy, p$K, pos = 3, cex = 0.8)
  invisible(x)
}
