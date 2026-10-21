# ==============================================================================
# Skewness and kurtosis tests of a fitted mixture
# ==============================================================================

#' Skewness and Kurtosis Tests of Model Fit
#'
#' @description
#' Checks whether a fitted mixture of normal distributions reproduces the
#' skewness and kurtosis of the data. A mixture of normals can be skewed and
#' heavy-tailed overall, through the spread of its classes; what it cannot
#' produce is a shape its classes do not add up to. A data set whose shape the
#' fitted model does not reproduce is evidence that the number of classes, or
#' the normality of the outcomes within each class, is wrong (Muthén, 2003).
#'
#' Each statistic is computed on the data and on `nsim` data sets simulated
#' from the fitted model; the simulated values give the statistic's
#' distribution under the model, and the p-value is two-sided: twice the
#' smaller of the shares of simulated values at or above, and at or below, the
#' sample value.
#'
#' The fitted parameters are estimates, not known values, and that uncertainty
#' belongs in the reference distribution. Each simulated data set is therefore
#' generated from its own parameter draw, taken from the normal approximation
#' to the estimates' sampling distribution (centred on the estimates, with the
#' inverse of the observed information as its covariance, on a scale where
#' variances are log standard deviations and class sizes log-ratios), with
#' classes drawn from the drawn class sizes and the estimation sample size.
#' Holding the parameters at the estimates instead understates the spread of
#' the skewness statistics markedly, because a mixture's skewness moves with
#' the separation and sizes of its classes, which are the least precisely
#' estimated parts of the model; the kurtosis statistics are dominated by
#' sampling noise and change little.
#'
#' @details
#' The statistics, for complete data:
#' \itemize{
#'   \item univariate: the skewness \eqn{m_3 / s^3} and the excess kurtosis
#'     \eqn{m_4 / s^4 - 3} of each outcome, \eqn{m_r} its \eqn{r}-th central
#'     sample moment and \eqn{s^2} its sample variance on \eqn{n - 1} degrees
#'     of freedom;
#'   \item bivariate and multivariate: Mardia's (1970) skewness
#'     \eqn{b_{1,p} = n^{-2} \sum_{i,j} (d_i^\top S^{-1} d_j)^3} and kurtosis
#'     \eqn{b_{2,p} = n^{-1} \sum_i (d_i^\top S^{-1} d_i)^2} of each pair of
#'     outcomes and of all of them, \eqn{d_i} the centred case and \eqn{S} the
#'     covariance matrix on \eqn{n} degrees of freedom. Under multivariate
#'     normality \eqn{b_{2,p}} is about \eqn{p(p + 2)}.
#' }
#' The p-values carry the simulation's own noise; with the default 100 data
#' sets they are good to about two decimals.
#'
#' @param fit A model with continuous outcomes from [`fit_gmm()`],
#'   [`fit_lcga()`] or [`fit_mixture()`], fitted to complete data.
#' @param nsim Number of simulated data sets.
#' @param seed `NULL`, or an integer for the simulation; a supplied seed
#'   leaves the caller's random number state as it was.
#'
#' @return A list of class `"skew_kurtosis_test"` with data frames
#'   `multivariate`, `univariate` and `bivariate`, each with columns
#'   `statistic` (`"skewness"` or `"kurtosis"`), the variable(s), `sample`,
#'   `mean`, `sd` and `p_value`.
#'
#' @references
#' Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with
#' applications. *Biometrika*, *57*(3), 519-530.
#'
#' Muthén, B. (2003). Statistical and substantive checking in growth mixture
#' modeling: Comment on Bauer and Curran (2003). *Psychological Methods*,
#' *8*(3), 369-377.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 300
#' cls <- rbinom(n, 1, 0.5)
#' y <- outer(ifelse(cls == 1, 3, 1), rep(1, 4)) +
#'   outer(ifelse(cls == 1, 1, 0.2), 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
#' fit <- fit_gmm(y, n_classes = 2, times = 4, random_effects = "none",
#'                n_init = 3, random_state = 1)
#' skew_kurtosis_test(fit, nsim = 20, seed = 1)
#' }
#' @export
skew_kurtosis_test <- function(fit, nsim = 100, seed = NULL) {
  if (!inherits(fit, "mixture_model") || is.null(fit$data))
    stop("`fit` must be a mixture model fitted by fit_gmm(), fit_lcga() or ",
         "fit_mixture().", call. = FALSE)
  Y <- as.matrix(fit$data)
  if (!is.numeric(Y) || anyNA(Y))
    stop("The test needs complete continuous outcomes; this fit's data are ",
         if (anyNA(Y)) "incomplete." else "not numeric.", call. = FALSE)
  fam <- class(fit$mm)[1L]
  if (!fam %in% c("structured_normal", "gaussian_diag", "gaussian_full",
                  "gaussian_unit", "gaussian_spherical", "gaussian_tied",
                  "gaussian_model"))
    stop("The test is for mixtures of normal outcomes; this fit's ",
         "measurement model is `", fam, "`.", call. = FALSE)
  nsim <- as.integer(nsim)
  if (is.na(nsim) || nsim < 2L)
    stop("`nsim` must be at least 2.", call. = FALSE)
  .simulate_refuse(fit)
  if (!is.null(fit$sm) && .supplies_class_probs(fit$sm))
    stop("The test is for a model without predictors of class membership; ",
         "fit the model without them first.", call. = FALSE)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  th <- .step1_pack(fit)
  if (is.null(th))
    stop("The test does not cover this model's parameterisation.",
         call. = FALSE)
  V <- .step1_variance(fit, Y, th, .step1_case_scores(fit, Y, th, rep(1, nrow(Y))),
                       rep(1, nrow(Y)), "hessian")$V
  # A semi-definite V (a boundary estimate) draws along its non-null directions.
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  R <- t(e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(th)))

  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      stats::runif(1)
    saved <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv))
    set.seed(seed)
  }
  obs <- .sk_stats(Y)
  n   <- nrow(Y); K <- fit$n_components
  S   <- vapply(seq_len(nsim), function(s) {
    ms  <- .step1_unpack(fit, th + drop(stats::rnorm(length(th)) %*% R))
    cls <- sample.int(K, n, replace = TRUE, prob = ms$weights)
    G   <- as.matrix(generate_synthetic_data(ms$mm, cls, n))
    colnames(G) <- colnames(Y)
    .sk_stats(G)$value
  }, numeric(nrow(obs)))
  S   <- matrix(S, nrow(obs))
  obs$mean <- rowMeans(S)
  obs$sd   <- apply(S, 1, stats::sd)
  obs$p_value <- vapply(seq_len(nrow(obs)), function(r)
    min(1, 2 * min(mean(S[r, ] >= obs$value[r]), mean(S[r, ] <= obs$value[r]))),
    numeric(1))
  names(obs)[names(obs) == "value"] <- "sample"
  split1 <- function(level) {
    d <- obs[obs$level == level, setdiff(names(obs), "level"), drop = FALSE]
    rownames(d) <- NULL
    d
  }
  out <- list(multivariate = split1("multivariate"),
              univariate = split1("univariate"),
              bivariate = split1("bivariate"), nsim = nsim)
  class(out) <- "skew_kurtosis_test"
  out
}

# Every statistic of one data matrix, in a fixed order: the multivariate pair,
# then each outcome's skewness and kurtosis, then each pair's.
.sk_stats <- function(Y) {
  p  <- ncol(Y); nm <- colnames(Y)
  mardia <- function(Z) {
    d <- sweep(Z, 2, colMeans(Z))
    S <- crossprod(d) / nrow(Z)
    D <- d %*% solve(S, t(d))
    c(mean(D^3), mean(diag(D)^2))
  }
  uni <- function(x) {
    d  <- x - mean(x)
    s2 <- sum(d^2) / (length(x) - 1)
    c(mean(d^3) / s2^1.5, mean(d^4) / s2^2 - 3)
  }
  rows <- list(data.frame(level = "multivariate",
                          statistic = c("skewness", "kurtosis"),
                          variables = "all", value = mardia(Y)))
  u <- vapply(seq_len(p), function(j) uni(Y[, j]), numeric(2))
  rows[[2]] <- data.frame(level = "univariate",
                          statistic = rep(c("skewness", "kurtosis"), each = p),
                          variables = rep(nm, 2), value = c(u[1, ], u[2, ]))
  if (p > 2L) {
    pr <- utils::combn(p, 2)
    b  <- vapply(seq_len(ncol(pr)), function(k) mardia(Y[, pr[, k]]), numeric(2))
    lab <- paste(nm[pr[1, ]], "and", nm[pr[2, ]])
    rows[[3]] <- data.frame(level = "bivariate",
                            statistic = rep(c("skewness", "kurtosis"), each = ncol(pr)),
                            variables = rep(lab, 2), value = c(b[1, ], b[2, ]))
  }
  do.call(rbind, rows)
}

#' @export
print.skew_kurtosis_test <- function(x, digits = 3, ...) {
  cat(sprintf("\nSkewness and kurtosis tests of model fit (%d simulated data sets)\n",
              x$nsim))
  for (lv in c("multivariate", "univariate", "bivariate")) {
    d <- x[[lv]]
    if (!nrow(d)) next
    cat(sprintf("\n%s%s\n", toupper(substr(lv, 1, 1)), substring(lv, 2)))
    d[c("sample", "mean", "sd")] <- lapply(d[c("sample", "mean", "sd")], round, digits)
    d$p_value <- round(d$p_value, 2)
    print(d, row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}
