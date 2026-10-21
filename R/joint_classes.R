# ==============================================================================
# Two latent class variables in one model
# ==============================================================================
#
# A loglinear model with two latent class variables (Hagenaars 1993): `ca`,
# measured by its own items, and `cs`, measured by its own, with the
# association written as a multinomial logit of `cs` on `ca`. That regression
# is exactly a two-occasion LTA's by-origin transition regression, so the fit
# rides fit_lta(): occasion 1 carries ca's items, occasion 2 cs's, each padded
# with all-NA columns to the longer list, measurement free by occasion, and a
# smaller `cs` written as a padded status whose destination column is
# forbidden. The padding is bookkeeping only and must not reach the user: its
# parameters touch no data, so they are taken out of the count (and every
# criterion built on it), and the report speaks of `ca` and `cs`, never of
# occasions or transitions.

#' Fit a model with two associated latent class variables
#'
#' Fits two latent class variables at once, each measured by its own items:
#' `ca` with `n_ca` classes and `cs` with `n_cs` classes, linked by a
#' multinomial logistic regression of `cs` on `ca` (Hagenaars, 1993). The
#' association is reported against the last class of each variable: the
#' intercept of `cs` class j is its log-odds against the last `cs` class among
#' members of the last `ca` class, and the coefficient of `ca` class k is how
#' much membership in `ca` class k changes that log-odds. The joint class
#' probabilities and the probabilities of each `cs` class given each `ca`
#' class come with it.
#'
#' Covariates can predict either variable. On `cs` each covariate acts in
#' common across the `ca` classes or separately within each, set per
#' covariate by `cs_effects`; the second form is a covariate whose effect on
#' `cs` is moderated by `ca`.
#'
#' @param ca_items,cs_items The items measuring each latent variable, one row
#'   per case, the same cases in the same order. Missing values are allowed.
#' @param n_ca,n_cs Integer. Number of classes of each variable. `n_ca` must be
#'   at least `n_cs`; when it is smaller, swap the roles of the two variables.
#' @param measurement_ca,measurement_cs `"binary"` (items coded 0/1) or
#'   `"ordinal"` (ordered categories coded from 0 or from 1). A binary and an
#'   ordinal side are fitted together as ordinal items, which for a two-category
#'   item is the same model.
#' @param predictors_ca A data frame of covariates predicting `ca` class
#'   membership, or `NULL`.
#' @param predictors_cs A data frame of covariates predicting `cs` class
#'   membership, or `NULL`.
#' @param cs_effects How each covariate in `predictors_cs` acts: `"common"`
#'   (one slope per `cs` class, the default), `"by_ca"` (a slope per `ca`
#'   class for every covariate), or a named vector such as
#'   `c(female = "by_ca")`, the unnamed covariates acting in common.
#' @param n_steps `1` (the default) fits everything at once. `2` is the
#'   two-step estimator (Bakk & Kuha, 2018): step 1 fits each variable's
#'   measurement model on its own, as two separate latent class analyses;
#'   step 2 holds those item parameters fixed and estimates the association
#'   and the covariate effects. The classes then mean the same thing in every
#'   structural model compared, and the standard errors carry step 1's
#'   uncertainty. The fit keeps step 1 as `step1`.
#' @param smoothing,bayes_constants,n_init,random_state,n_cores,standard_errors
#'   As in [`fit_lta()`].
#' @param ... Further arguments passed to [`fit_lta()`], such as `weights`.
#'
#' @return An object of class `joint_classes` (inheriting from `lta_model`, so
#'   [`lr_test()`] and the information criteria apply), with `n_ca`, `n_cs`,
#'   `association` (from [`joint_coefficients()`]), `joint` (the `n_ca x n_cs`
#'   joint class probabilities) and `conditional` (P(cs | ca), rows summing to
#'   one).
#'
#' @references
#' Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between latent
#' classes and external variables. *Psychometrika*, 83(4), 871-892.
#' \doi{10.1007/s11336-017-9592-7}
#'
#' Hagenaars, J. A. (1993). *Loglinear Models with Latent Variables*. Sage.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 600
#' ca <- sample(1:2, n, TRUE)
#' cs <- ifelse(runif(n) < ifelse(ca == 1, 0.8, 0.2), 1L, 2L)
#' A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(ca == 1, 0.85, 0.15)))
#' S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cs == 1, 0.85, 0.15)))
#' fit <- fit_joint_classes(A, S, n_ca = 2, n_cs = 2, n_init = 5,
#'                          random_state = 1)
#' fit
#' }
#' @export
fit_joint_classes <- function(ca_items, cs_items,
                              n_ca = 2, n_cs = 2,
                              measurement_ca = c("binary", "ordinal"),
                              measurement_cs = c("binary", "ordinal"),
                              predictors_ca = NULL,
                              predictors_cs = NULL,
                              cs_effects = "common",
                              n_steps = 1,
                              smoothing = 1.0,
                              bayes_constants = NULL,
                              n_init = 20,
                              random_state = NULL,
                              n_cores = .default_n_cores(),
                              standard_errors = TRUE,
                              ...) {
  measurement_ca <- match.arg(measurement_ca)
  measurement_cs <- match.arg(measurement_cs)
  n_steps <- as.integer(n_steps)[1L]
  if (is.na(n_steps) || !(n_steps %in% 1:2))
    stop("`n_steps` must be 1 or 2.", call. = FALSE)
  A <- as.matrix(ca_items); S <- as.matrix(cs_items)
  if (nrow(A) != nrow(S))
    stop("`ca_items` and `cs_items` must have the same number of rows.",
         call. = FALSE)
  n_ca <- as.integer(n_ca); n_cs <- as.integer(n_cs)
  if (n_ca < 2L || n_cs < 2L)
    stop("`n_ca` and `n_cs` must each be at least 2.", call. = FALSE)
  # A smaller `ca` would need a padded class that nobody starts in, and the
  # initial-status probabilities have no zero to hold it at.
  if (n_ca < n_cs)
    stop("`n_ca` must be at least `n_cs`; swap the roles of the two ",
         "variables.", call. = FALSE)

  ordinal <- "ordinal" %in% c(measurement_ca, measurement_cs)
  if (ordinal) {
    # The ordinal emission reads categories from 1; a 0/1 item becomes 1/2,
    # the same two-category model.
    up <- function(M) apply(M, 2, function(x)
      if (isTRUE(min(x, na.rm = TRUE) == 0)) x + 1 else x)
    A <- up(A); S <- up(S)
  }
  Ja <- ncol(A); Js <- ncol(S); J <- max(Ja, Js); K <- n_ca
  pad <- function(M) cbind(M, matrix(NA_real_, nrow(M), J - ncol(M)))
  X <- cbind(pad(A), pad(S))

  forb <- NULL
  if (n_cs < K) { forb <- matrix(FALSE, K, K); forb[, (n_cs + 1L):K] <- TRUE }

  effects <- "common"
  if (!is.null(predictors_cs)) {
    predictors_cs <- as.data.frame(predictors_cs)
    effects <- .joint_cs_effects(cs_effects, names(predictors_cs))
  }

  # The two-step's step 1 is the two separate LCAs: with no association and
  # no covariates the joint likelihood factorises into them, which is the
  # transition-free fit (`.twostep_independent`, R/lta.R).
  fit <- fit_lta(X, n_statuses = K, times = 2,
                 n_steps = n_steps, .twostep_independent = n_steps == 2L,
                 measurement = if (ordinal) "ordinal" else "binary",
                 measurement_invariance = "none",
                 forbidden_transitions = forb,
                 predictors_initial = predictors_ca,
                 predictors_transition = predictors_cs,
                 transition_effects = effects,
                 smoothing = smoothing, bayes_constants = bayes_constants,
                 n_init = n_init, random_state = random_state,
                 n_cores = n_cores, standard_errors = standard_errors, ...)

  fit$n_ca <- n_ca; fit$n_cs <- n_cs
  fit$n_items_ca <- Ja; fit$n_items_cs <- Js
  n_pad <- .joint_padding_n_params(fit)
  fit$n_params <- fit$n_params - n_pad
  fit$metrics$n_params <- fit$n_params
  fresh <- .lta_metrics(fit)
  for (nm in intersect(c("aic", "bic", "caic", "awe", "aic3", "icl", "sabic"),
                       names(fresh)))
    fit$metrics[[nm]] <- fresh[[nm]]

  P <- fit$delta * fit$tau[[1]][, seq_len(n_cs), drop = FALSE]
  dimnames(P) <- list(paste0("CA#", seq_len(n_ca)), paste0("CS#", seq_len(n_cs)))
  fit$joint <- P
  fit$conditional <- P / rowSums(P)
  class(fit) <- c("joint_classes", class(fit))
  fit$association <- joint_coefficients(fit)
  fit
}

# `cs_effects` in fit_lta()'s terms: "by_ca" for every covariate is the named
# vector naming them all, which is the same model as a separate regression per
# origin (the origin dummies already give each origin its own intercept).
.joint_cs_effects <- function(cs_effects, covs) {
  if (is.null(names(cs_effects))) {
    eff <- match.arg(cs_effects, c("common", "by_ca"))
    if (eff == "common") return("common")
    return(stats::setNames(rep("by_origin", length(covs)), covs))
  }
  bad <- setdiff(unname(cs_effects), c("common", "by_ca"))
  if (length(bad))
    stop("`cs_effects` values must be \"common\" or \"by_ca\".", call. = FALSE)
  unknown <- setdiff(names(cs_effects), covs)
  if (length(unknown))
    stop("`cs_effects` names covariates not in `predictors_cs`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  stats::setNames(ifelse(cs_effects == "by_ca", "by_origin", "common"),
                  names(cs_effects))
}

# Free item parameters the padding added: the padded item positions on each
# side, and every item of the padded `cs` classes. Each cell costs what its
# column costs per class (one for a binary item, categories - 1 for an ordinal
# one). The padded destination column's transition cells are forbidden, so
# fit_lta() never counted them; with covariates on `cs` the regression has no
# category for it (.lta_tau_keep()), so the same holds.
.joint_padding_n_params <- function(fit) {
  K  <- fit$n_statuses
  Ja <- fit$n_items_ca; Js <- fit$n_items_cs; J <- max(Ja, Js)
  cats <- fit$mm$models[[1]]$cats
  per_item <- if (is.null(cats)) rep(1L, J) else as.integer(cats) - 1L
  pad_a <- if (Ja < J) K * sum(per_item[(Ja + 1L):J]) else 0L
  pad_s <- if (Js < J) K * sum(per_item[(Js + 1L):J]) else 0L
  pad_k <- (K - fit$n_cs) * sum(per_item[seq_len(Js)])
  pad_a + pad_s + pad_k
}

#' Association coefficients of a two-latent-variable model
#'
#' The multinomial logit of `cs` on `ca` from [`fit_joint_classes()`], in the
#' reference-class form: every log-odds is against the last `cs` class, and
#' `ca` enters as indicators against the last `ca` class. Rows are the
#' intercept (`term = "Intercept"`, the log-odds among the last `ca` class),
#' the `ca` effects (`term = "CA#k"`), and each covariate on `cs`: one row for
#' a common slope, or, for a covariate acting by `ca` class, its slope in the
#' last `ca` class (`ca = NA`) and the difference from it in each other class.
#'
#' @param fit A `joint_classes` fit.
#' @return A data frame with columns `cs`, `term`, `ca`, `estimate`, `se`,
#'   `z`, `p`.
#' @export
joint_coefficients <- function(fit) {
  stopifnot(inherits(fit, "joint_classes"))
  K1 <- fit$n_ca; K2 <- fit$n_cs
  V  <- fit$se$vcov
  np <- if (is.null(V)) 0L else ncol(V)
  eta <- .joint_eta(fit)
  cov <- !is.null(fit$tau_beta)
  cn  <- if (cov) colnames(.lta_tau_design(fit, 1L)) else character()

  # log-odds of cs class j against the last, for origin a
  lo <- function(a, j, term) {
    e1 <- eta(a, j, term); e2 <- eta(a, K2, term)
    list(v = e1$v - e2$v, r = e1$r - e2$r)
  }
  minus <- function(x, y) list(v = x$v - y$v, r = x$r - y$r)

  rows <- list()
  add <- function(j, term, ca, x)
    rows[[length(rows) + 1L]] <<- data.frame(
      cs = sprintf("CS#%d", j), term = term, ca = ca, estimate = x$v,
      se = if (np) sqrt(max(drop(crossprod(x$r, V %*% x$r)), 0)) else NA_real_)
  covs <- if (cov) cn[!grepl("^(Intercept$|from:)", cn)] else character()
  for (j in seq_len(K2 - 1L)) {
    add(j, "Intercept", NA_character_, lo(K1, j, "Intercept"))
    for (k in seq_len(K1 - 1L))
      add(j, sprintf("CA#%d", k), sprintf("CA#%d", k),
          minus(lo(k, j, "Intercept"), lo(K1, j, "Intercept")))
    for (v in covs) {
      add(j, v, NA_character_, lo(K1, j, v))
      if (cov && any(sprintf("from:1 x %s", v) == cn))
        for (k in seq_len(K1 - 1L))
          add(j, v, sprintf("CA#%d", k), minus(lo(k, j, v), lo(K1, j, v)))
    }
  }
  out <- do.call(rbind, rows)
  out$z <- out$estimate / out$se
  out$p <- 2 * stats::pnorm(-abs(out$z))
  out
}

# eta(a, b, term): the linear predictor of destination b for origin a, as a
# value and a row over the packed parameter vector, so every reported
# coefficient is a difference of etas and its SE one quadratic form. Without
# covariates each origin's row is a softmax anchored at its last allowed
# destination; the block's `ref` lists the free destinations in packing order.
# With covariates the regression's rows are its own destinations
# (.lta_tau_keep(): the `cs` classes, the padded ones dropped), the last the
# reference.
.joint_eta <- function(fit) {
  K  <- fit$n_statuses
  V  <- fit$se$vcov
  np <- if (is.null(V)) 0L else ncol(V)
  if (!is.null(fit$tau_beta)) {
    B  <- fit$tau_beta[[1]]; cn <- colnames(.lta_tau_design(fit, 1L))
    blk <- Find(function(b) grepl("^tau_beta", b$name), fit$se$blocks)
    keep <- .lta_tau_keep(fit)
    return(function(a, b, term) {
      nm <- if (term == "Intercept") c("Intercept", if (a < K) sprintf("from:%d", a))
            else c(term, intersect(sprintf("from:%d x %s", a, term), cn))
      r <- numeric(np)
      d <- match(b, keep)
      if (d == length(keep)) return(list(v = 0, r = r))
      idx <- match(nm, cn)
      if (!is.null(blk)) r[blk$cols[(d - 1L) * length(cn) + idx]] <- 1
      list(v = sum(B[d, idx]), r = r)
    })
  }
  lt <- log(pmax(fit$tau[[1]], 1e-300))
  function(a, b, term) {
    r <- numeric(np)
    blk <- Find(function(x) x$name == sprintf("tau(1)[from %d]", a),
                fit$se$blocks)
    if (!is.null(blk) && b %in% blk$ref) r[blk$cols[match(b, blk$ref)]] <- 1
    list(v = lt[a, b], r = r)
  }
}

#' @export
print.joint_classes <- function(x, digits = 3, ...) {
  cat(sprintf("Two latent class variables: ca (%d classes, %d items), cs (%d classes, %d items)\n",
              x$n_ca, x$n_items_ca, x$n_cs, x$n_items_cs))
  cat(sprintf("Log-likelihood %.3f, %d free parameters, BIC %.3f\n\n",
              x$metrics$ll, as.integer(x$n_params), x$metrics$bic))
  if (identical(x$n_steps, 2L))
    cat("Two-step: item parameters fixed at the separate latent class analyses\n\n")
  cat("Association: log-odds of cs against its last class\n")
  a <- x$association
  a$ca[is.na(a$ca)] <- ""
  print(format(a, digits = digits), row.names = FALSE)
  cat("\nJoint class probabilities P(ca, cs)\n")
  print(round(x$joint, digits))
  cat("\nConditional probabilities P(cs | ca)\n")
  print(round(x$conditional, digits))
  invisible(x)
}

# The inherited LTA summary speaks of occasions and transitions; this model has
# neither, so its summary is its own print.
#' @export
summary.joint_classes <- function(object, digits = 3, ...) {
  print.joint_classes(object, digits = digits, ...)
}
