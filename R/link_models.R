# ------------------------------------------------------------------------------
# Linking separately fitted mixture models by the bias-adjusted three-step.
# ------------------------------------------------------------------------------
#
# The question this answers: a latent class model at one occasion and a growth
# mixture model at a later one (or any two mixture models, of any family and
# any number of classes), and does membership in the first predict membership
# in the second? Each model is step 1 on its own. Step 2 classifies every case
# under each model and forms that model's classification-error table. Step 3
# regresses the later latent class on the earlier one, with the error tables
# held fixed, exactly as the three-step fit_lta() does for one model repeated
# over occasions.
#
# WHY THIS IS NOT fit_lta(n_steps = 3). That function's step 3 is itself a
# fit_lta() call, and a latent transition model has one number of statuses at
# every occasion: its forward-backward recursion, its packing and its
# transition regressions are all built on a K x K matrix. Here the matrices are
# K_t x K_{t+1}. Padding the smaller side with statuses nobody can enter would
# need structural zeros in the transitions, which a transition regression on
# covariates cannot hold. Step 3 has only the structural parameters left once
# the error tables are fixed, a dozen or so, so it is maximised directly.
#
# THE PARAMETERISATION, last class the reference throughout. The first
# occasion's class is a multinomial logit on `predictors_initial`. Each later
# occasion's class is a multinomial logit on dummies for the previous
# occasion's class (last class the reference) and on `predictors_transition`,
# one slope per destination shared by every origin. With no covariates and two
# occasions that is the saturated model for the K_1 x K_2 table.

#' Link separately fitted mixture models by the three-step
#'
#' Relates the latent classes of two or more mixture models fitted on the same
#' cases -- for example a latent class analysis of early attitudes and a growth
#' mixture model of later achievement -- with the bias-adjusted three-step
#' estimator. Each fitted model is step 1. Step 2 assigns every case to its
#' most likely class under each model and tabulates how often that assignment
#' misses the true class. Step 3 fits a multinomial logistic regression of each
#' occasion's latent class on the previous occasion's, holding those error
#' tables fixed, so the classification error does not attenuate the
#' association (Vermunt, 2010; Bakk, Tekle & Vermunt, 2013).
#'
#' The models may differ in family and in number of classes. Occasion 1's
#' class is regressed on `predictors_initial`; each later occasion's class on
#' the previous occasion's class and on `predictors_transition`, with one slope
#' per destination class shared by every origin class. The last class is the
#' reference throughout.
#'
#' @param models A list of two or more fitted models (from [fit_mixture()],
#'   [fit_gmm()] or similar), in occasion order, each fitted on the same cases
#'   in the same row order. Names, if given, label the occasions.
#' @param predictors_initial,predictors_transition Optional covariates for the
#'   first occasion's class and for the later occasions' classes: a numeric
#'   vector, matrix or data frame with one row per case.
#' @param correction `"ML"` (default) holds the classification error fixed in
#'   step 3; `"none"` treats the assigned classes as the true ones, the naive
#'   classify-analyse baseline.
#' @param weights Optional case weights, treated as frequencies.
#' @param n_init Number of starting points for step 3; the first is all
#'   coefficients zero, the rest random.
#' @param random_state Optional seed for the random starts.
#' @param standard_errors Logical; compute standard errors from the observed
#'   information of step 3. They treat the error tables as known, so they do
#'   not include the uncertainty of the step-1 estimates.
#'
#' @return An object of class `linked_model` with the step-3 log-likelihood
#'   (`loglik`), the number of parameters, BIC, the coefficients and their
#'   covariance matrix, the implied initial class proportions (`initial`) and
#'   transition probabilities (`transitions`, one K_t x K_{t+1} matrix per pair
#'   of occasions, averaged over cases), the classification-error tables (rows
#'   the assigned class, columns the true class) and the modal assignments.
#'
#' @references
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two improved
#' three-step approaches. *Political Analysis*, 18(4), 450-469.
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the association
#' between latent class membership and external variables using bias-adjusted
#' three-step approaches. *Sociological Methodology*, 43(1), 272-311.
#'
#' Nylund-Gibson, K., Grimm, R., Quirk, M., & Furlong, M. (2014). A latent
#' transition mixture model using the three-step specification. *Structural
#' Equation Modeling*, 21(3), 439-454.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 400
#' s1 <- sample(1:3, n, replace = TRUE)
#' s2 <- ifelse(runif(n) < c(0.8, 0.5, 0.2)[s1], 1, 2)
#' X1 <- sapply(1:5, function(j) rbinom(n, 1, c(0.9, 0.5, 0.1)[s1]))
#' X2 <- sapply(1:4, function(j) rnorm(n, c(0, 3)[s2]))
#' a <- fit_mixture(X1, n_classes = 3, measurement = "binary", n_init = 5)
#' b <- fit_mixture(X2, n_classes = 2, measurement = "continuous", n_init = 5)
#' link <- link_models(list(early = a, late = b))
#' link
#' }
#' @export
link_models <- function(models,
                        predictors_initial    = NULL,
                        predictors_transition = NULL,
                        correction            = c("ML", "none"),
                        weights               = NULL,
                        n_init                = 10,
                        random_state          = NULL,
                        standard_errors       = TRUE) {
  correction <- match.arg(correction)
  if (!is.list(models) || inherits(models, "mixture_model") || length(models) < 2L)
    stop("`models` must be a list of at least two fitted models, one per ",
         "occasion.", call. = FALSE)
  if (any(vapply(models, inherits, TRUE, "lta_model")))
    stop("`models` must be cross-sectional mixture models; a latent ",
         "transition model is linked over its occasions by ",
         "`fit_lta(n_steps = 3)`.", call. = FALSE)

  post <- lapply(models, class_assignments, type = "posterior")
  n <- nrow(post[[1L]])
  if (any(vapply(post, nrow, 1L) != n))
    stop("The models were fitted on different numbers of cases. Fit every ",
         "model on the same cases, in the same row order.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, n)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0))
    stop("`weights` must be one non-negative number per case.", call. = FALSE)

  labels <- names(models)
  if (is.null(labels) || any(!nzchar(labels)))
    labels <- paste("Time", seq_along(models))

  # Step 2, the rule the three-step fit_lta() uses: modal assignment, one
  # table per model, rows the assigned class and columns the true one.
  W <- vapply(post, max.col, integer(n), ties.method = "first")
  D <- lapply(post, function(p) {
    d <- .classification_error(p, assignment = "modal", weights = weights)$D
    if (correction == "none") d[] <- diag(ncol(p))
    d
  })

  prep <- function(x, label) {
    if (is.null(x)) return(NULL)
    .lta_design(x, n, label)[, -1L, drop = FALSE]
  }
  fit <- .link_step3(W, D, prep(predictors_initial, "predictors_initial"),
                     prep(predictors_transition, "predictors_transition"),
                     weights, n_init = n_init, random_state = random_state,
                     standard_errors = standard_errors, labels = labels)
  fit$correction <- correction
  fit$models     <- models
  fit
}

# Step 3 on the assigned classes. `W` is n x T (modal classes), `D` a list of
# the T error tables (rows assigned, columns true), `Zd` / `Zt` the covariate
# matrices without an intercept, or NULL. Separate from link_models() so that
# a step 1 and step 2 obtained elsewhere can be handed in directly.
.link_step3 <- function(W, D, Zd = NULL, Zt = NULL, weights = NULL,
                        n_init = 10, random_state = NULL,
                        standard_errors = TRUE, labels = NULL) {
  W  <- as.matrix(W)
  n  <- nrow(W)
  Tn <- ncol(W)
  K  <- vapply(D, ncol, 1L)
  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(labels)) labels <- paste("Time", seq_len(Tn))
  nd <- if (is.null(Zd)) 0L else ncol(Zd)
  nt <- if (is.null(Zt)) 0L else ncol(Zt)
  cov_names <- function(Z, nz) if (nz)
    colnames(Z) %||% paste0("x", seq_len(nz)) else character(0)

  # The likelihood reads a case only through its assigned classes and its
  # covariate values, so cases that agree on both are one row with the summed
  # weight: exact, and a few dozen rows instead of the sample.
  key  <- do.call(paste, c(as.data.frame(cbind(W, Zd, Zt)), sep = "\r"))
  grp  <- match(key, unique(key))
  keep <- match(seq_len(max(grp)), grp)
  w    <- as.vector(rowsum(weights, grp))
  Wp   <- W[keep, , drop = FALSE]
  P    <- length(w)
  Zdp  <- cbind(matrix(1, P, 1L), if (nd) Zd[keep, , drop = FALSE])
  Zt_only <- if (nt) Zt[keep, , drop = FALSE]
  logD <- lapply(seq_len(Tn), function(t)
    log(pmax(D[[t]], 1e-300))[Wp[, t], , drop = FALSE])

  # The parameter vector: the initial block, (K_1 - 1) x (1 + nd), then per
  # transition a (K_{t+1} - 1) x (1 + (K_t - 1) + nt) block of intercept,
  # origin dummies and slopes, each block stored column by column.
  dims <- c(list(c(K[1L] - 1L, 1L + nd)),
            lapply(seq_len(Tn - 1L), function(t)
              c(K[t + 1L] - 1L, 1L + (K[t] - 1L) + nt)))
  len  <- vapply(dims, prod, 1)
  end  <- cumsum(len)
  block <- function(par, b) matrix(par[(end[b] - len[b] + 1):end[b]], dims[[b]][1L])

  log_softmax0 <- function(eta) {
    eta <- cbind(eta, 0)
    eta - logsumexp(eta)
  }
  # Log initial probabilities, P x K_1; log transitions out of origin k at
  # transition t, P x K_{t+1}.
  log_init <- function(par) log_softmax0(Zdp %*% t(block(par, 1L)))
  log_tau  <- function(par, t, k) {
    B   <- block(par, t + 1L)
    Kt  <- K[t]
    eta <- matrix(B[, 1L], P, nrow(B), byrow = TRUE)
    if (k < Kt) eta <- eta + matrix(B[, 1L + k], P, nrow(B), byrow = TRUE)
    if (nt) eta <- eta + Zt_only %*% t(B[, Kt + seq_len(nt), drop = FALSE])
    log_softmax0(eta)
  }

  # Forward recursion over the occasions: la[i, s] is the log of the joint
  # probability of the assignments so far and true class s now.
  ll_pattern <- function(par) {
    la <- log_init(par) + logD[[1L]]
    for (t in seq_len(Tn - 1L)) {
      nxt <- matrix(-Inf, P, K[t + 1L])
      terms <- lapply(seq_len(K[t]), function(k) la[, k] + log_tau(par, t, k))
      for (j in seq_len(K[t + 1L]))
        nxt[, j] <- logsumexp(vapply(terms, function(x) x[, j], numeric(P)))
      la <- nxt + logD[[t + 1L]]
    }
    logsumexp(la)
  }
  ll <- function(par) {
    v <- sum(w * ll_pattern(par))
    if (is.finite(v)) v else -1e300
  }

  npar <- sum(len)
  if (!is.null(random_state)) set.seed(random_state)
  starts <- c(list(numeric(npar)),
              lapply(seq_len(max(0L, n_init - 1L)), function(i) stats::rnorm(npar)))
  runs <- lapply(starts, function(s)
    stats::optim(s, ll, method = "BFGS",
                 control = list(fnscale = -1, maxit = 1000, reltol = 1e-14)))
  lls  <- vapply(runs, `[[`, 0, "value")
  best <- runs[[which.max(lls)]]$par

  # A few Newton steps on central differences take BFGS's stopping point to
  # the maximum; accepted only while they raise the log-likelihood.
  H <- NULL
  for (it in 1:5) {
    g <- .link_grad(ll, best)
    H <- .link_hessian(ll, best)
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    cand <- best - step
    if (ll(cand) < ll(best)) break
    best <- cand
    if (max(abs(step)) < 1e-10) break
  }
  loglik <- ll(best)

  # Names in the model's own terms: the destination class, then what it is
  # regressed on.
  cls <- function(t, k) sprintf("%s Class %d", labels[t], k)
  nm <- c(as.vector(outer(seq_len(K[1L] - 1L), c("(Intercept)", cov_names(Zd, nd)),
                          function(k, v) paste(cls(1L, k), "~", v))))
  for (t in seq_len(Tn - 1L)) {
    rhs <- c("(Intercept)", if (K[t] > 1L) cls(t, seq_len(K[t] - 1L)),
             cov_names(Zt, nt))
    nm <- c(nm, as.vector(outer(seq_len(K[t + 1L] - 1L), rhs,
                                function(k, v) paste(cls(t + 1L, k), "~", v))))
  }
  names(best) <- nm

  vc <- NULL
  if (standard_errors) {
    if (is.null(H)) H <- .link_hessian(ll, best)
    vc <- tryCatch(solve(-H), error = function(e) NULL)
    if (!is.null(vc)) dimnames(vc) <- list(nm, nm)
  }

  # Average implied proportions and transition probabilities over the cases.
  avg <- function(M) colSums(M * w) / sum(w)
  initial <- avg(exp(log_init(best)))
  names(initial) <- cls(1L, seq_len(K[1L]))
  transitions <- lapply(seq_len(Tn - 1L), function(t) {
    M <- t(vapply(seq_len(K[t]), function(k) avg(exp(log_tau(best, t, k))),
                  numeric(K[t + 1L])))
    dimnames(M) <- list(cls(t, seq_len(K[t])), cls(t + 1L, seq_len(K[t + 1L])))
    M
  })
  names(transitions) <- paste(labels[-Tn], "->", labels[-1L])

  structure(list(
    loglik               = loglik,
    n_parameters         = npar,
    n_obs                = sum(weights),
    bic                  = -2 * loglik + npar * log(sum(weights)),
    coefficients         = best,
    vcov                 = vc,
    se                   = if (is.null(vc)) NULL else sqrt(pmax(diag(vc), 0)),
    initial              = initial,
    transitions          = transitions,
    classification_error = D,
    modal                = W,
    n_classes            = K,
    labels               = labels,
    start_loglik         = lls),
    class = "linked_model")
}

.link_grad <- function(f, par, h = 1e-5) {
  vapply(seq_along(par), function(j) {
    e <- replace(numeric(length(par)), j, h)
    (f(par + e) - f(par - e)) / (2 * h)
  }, 0)
}

.link_hessian <- function(f, par, h = 1e-4) {
  p <- length(par)
  H <- matrix(0, p, p)
  f0 <- f(par)
  for (i in seq_len(p)) {
    ei <- replace(numeric(p), i, h)
    H[i, i] <- (f(par + ei) - 2 * f0 + f(par - ei)) / h^2
    for (j in seq_len(i - 1L)) {
      ej <- replace(numeric(p), j, h)
      H[i, j] <- H[j, i] <- (f(par + ei + ej) - f(par + ei - ej) -
                               f(par - ei + ej) + f(par - ei - ej)) / (4 * h^2)
    }
  }
  H
}

#' @export
print.linked_model <- function(x, digits = 3, ...) {
  cat(sprintf("Three-step link of %d mixture models (%s); classes %s\n",
              length(x$n_classes), paste(x$labels, collapse = ", "),
              paste(x$n_classes, collapse = " -> ")))
  cat(sprintf("Log-likelihood %.3f, %d parameters, BIC %.3f\n\n",
              x$loglik, x$n_parameters, x$bic))
  tab <- cbind(Estimate = x$coefficients,
               SE = if (is.null(x$se)) NA_real_ else x$se)
  print(round(tab, digits))
  cat("\nTransition probabilities (averaged over cases):\n")
  for (nm in names(x$transitions)) {
    cat(nm, "\n")
    print(round(x$transitions[[nm]], digits))
  }
  cat("\nStandard errors treat the classification-error tables as known.\n")
  invisible(x)
}

#' @export
coef.linked_model <- function(object, ...) object$coefficients

#' @export
vcov.linked_model <- function(object, ...) object$vcov

#' @export
logLik.linked_model <- function(object, ...)
  structure(object$loglik, df = object$n_parameters, nobs = object$n_obs,
            class = "logLik")
