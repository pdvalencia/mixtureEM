# ==============================================================================
# Wald test of linear constraints on named parameters
# ==============================================================================

#' Wald Test of Linear Constraints on Model Parameters
#'
#' @description
#' Tests a set of linear equality constraints on the estimated parameters of a
#' fitted model -- for example, that a growth-factor mean is the same in three
#' classes -- without refitting the model under the constraints. The
#' constraints name parameters by label; the labels are those of
#' [`growth_estimates()`] for a growth model and those of
#' [`vcov()`][vcov.mixture_model] for a covariate model (see Details).
#'
#' @details
#' Writing the constraints as \eqn{R\theta = r}, the statistic is
#' \deqn{W = (R\hat\theta - r)^\top (R \hat V R^\top)^{-1} (R\hat\theta - r),}
#' referred to a chi-squared distribution with as many degrees of freedom as
#' there are constraints (Wald, 1943). \eqn{\hat V} is the estimated covariance
#' matrix of the parameters.
#'
#' **Labels.** For a [`fit_gmm()`] or [`fit_lcga()`] fit, every row of
#' `growth_estimates(fit)` is one parameter, labelled
#' `"<parameter>[<item>, <class>]"`: `"growth_mean[intercept, 1]"`,
#' `"growth_variance[linear, 2]"`, `"residual_variance[t1, 1]"`, and so on.
#' A parameter constrained equal across classes appears once per class with
#' the same value, and equating two of its copies is not a constraint at all.
#' For a mixture model with a covariate structural model the labels are the
#' names of `vcov(fit)`, `"Class k:predictor"`, which are the multinomial-logit
#' coefficients relative to the class reported in its `ref_class` attribute.
#'
#' **Constraints.** Either a list, each element one constraint:
#' * a pair of labels, `c("a", "b")`, meaning `a - b = 0`;
#' * a named numeric vector of weights, `c(a = 1, b = -2)`, meaning
#'   `a - 2 b = rhs`;
#'
#' or a numeric matrix with one row per constraint whose column names are
#' labels. Three classes are equal on a parameter when two pairwise equalities
#' hold, so `list(c("a1", "a2"), c("a2", "a3"))` has 2 degrees of freedom; a
#' third equality, `c("a1", "a3")`, would add nothing and is refused as
#' redundant.
#'
#' **Which test to use.** `wald_test()` answers "are these parameters
#' equal (or equal to these values)?" for any parameters of the model.
#' [`analytical_wald_test()`] answers the narrower and common question "does
#' this covariate affect class membership at all?", and needs only the
#' covariate's name; on a covariate model it is the same statistic as
#' `wald_test()` with every one of that covariate's coefficients set to zero.
#' [`wald_omnibus_test()`] answers the same question from bootstrapped
#' standard errors, for samples too small to trust the analytical ones.
#'
#' Like any Wald test, this one assumes the estimates are near-normal, which
#' fails near the edge of the parameter space: a variance close to zero, or a
#' class with very few cases. A likelihood-ratio test of the constrained model
#' against the free one settles those cases.
#'
#' @param fit A model from [`fit_gmm()`] or [`fit_lcga()`], or a
#'   `mixture_model` with a covariate structural model.
#' @param constraints A list or matrix of constraints. See Details.
#' @param rhs Numeric, the right-hand side \eqn{r}: one value per constraint,
#'   or a single value recycled. Default `0`.
#' @param se For growth models, which covariance estimator to use: `"robust"`
#'   (default), `"hessian"` or `"opg"`, as in [`growth_estimates()`]. Ignored
#'   for covariate models, whose covariance is the one the model was fitted
#'   with (see [`vcov()`][vcov.mixture_model]).
#'
#' @return A one-row data frame of class `mixture_wald` with columns
#'   `Wald_Chi2`, `df`, `p_value` and `Method`, the name of the covariance
#'   estimator behind the statistic.
#'
#' @references
#' Wald, A. (1943). Tests of statistical hypotheses concerning several
#' parameters when the number of observations is large. *Transactions of the
#' American Mathematical Society*, *54*(3), 426-482.
#' \doi{10.1090/S0002-9947-1943-0012401-3}
#'
#' @seealso [`growth_estimates()`] for the labels and estimates of a growth
#'   model, [`analytical_wald_test()`], [`wald_omnibus_test()`].
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 400
#' cls <- rbinom(n, 1, 0.5)
#' icept <- ifelse(cls == 1, 3, 1) + rnorm(n, 0, 1)
#' slope <- ifelse(cls == 1, 1.0, 0.2) + rnorm(n, 0, 0.4)
#' y <- outer(icept, rep(1, 4)) + outer(slope, 0:3) + matrix(rnorm(n * 4, 0, 0.7), n, 4)
#' fit <- fit_gmm(y, times = 4, n_classes = 2, n_init = 10, random_state = 1)
#' # Do the two classes share their mean slope?
#' wald_test(fit, list(c("growth_mean[linear, 1]", "growth_mean[linear, 2]")))
#' }
#' @export
wald_test <- function(fit, constraints, rhs = 0,
                      se = c("robust", "hessian", "opg")) {
  se <- match.arg(se)
  pe <- .wald_parameters(fit, se)
  R  <- .wald_constraint_matrix(constraints, names(pe$estimate))
  q  <- nrow(R)
  if (!length(rhs) %in% c(1L, q))
    stop(sprintf("`rhs` must have length 1 or %d (one per constraint).", q),
         call. = FALSE)
  r <- rep_len(as.numeric(rhs), q)

  if (qr(R)$rank < q)
    stop(paste0("The constraints are redundant: at least one follows from the ",
                "others. Drop it; k classes are made equal on a parameter by ",
                "k - 1 equalities."), call. = FALSE)
  V <- pe$vcov
  if (anyNA(V[, colSums(R != 0) > 0]))
    stop("A constrained parameter has no standard error (see growth_estimates()).",
         call. = FALSE)
  d  <- as.vector(R %*% pe$estimate) - r
  RV <- R %*% V %*% t(R)
  W  <- tryCatch(as.numeric(crossprod(d, solve(RV, d))),
                 error = function(e) stop(paste0(
                   "The constrained combinations have a singular covariance, ",
                   "so the statistic is not defined. A constraint on a fixed ",
                   "parameter, or on two copies of one parameter shared by ",
                   "classes, does this."), call. = FALSE))

  out <- data.frame(Wald_Chi2 = W, df = q,
                    p_value = stats::pchisq(W, q, lower.tail = FALSE),
                    Method = pe$method, stringsAsFactors = FALSE)
  class(out) <- c("mixture_wald", "data.frame")
  out
}

# Named estimates and their covariance for whichever kind of fit this is.
.wald_parameters <- function(fit, se) {
  if (inherits(fit, c("gmm", "lcga"))) {
    ge <- growth_estimates(fit, se = se)
    V  <- attr(ge, "vcov")
    if (is.null(V))
      stop("This growth model has no covariance matrix (see growth_estimates()).",
           call. = FALSE)
    lab <- sprintf("%s[%s, %d]", ge$parameter, ge$item, as.integer(ge$class))
    est <- stats::setNames(ge$estimate, lab)
    dimnames(V) <- list(lab, lab)
    method <- c(robust = "Sandwich (robust)", hessian = "Observed information",
                opg = "Outer product of scores")[[se]]
    return(list(estimate = est, vcov = V, method = method))
  }
  if (!is.null(fit$sm) && inherits(fit$sm, "covariate")) {
    V   <- vcov(fit)
    ref <- attr(V, "ref_class")
    b   <- coef(fit, exponentiate = FALSE, ref_class = ref)
    D   <- ncol(b)
    free <- setdiff(seq_len(nrow(b)), ref)
    est <- as.vector(t(b[free, , drop = FALSE]))
    names(est) <- rownames(V)
    return(list(estimate = est, vcov = V,
                method = attr(V, "method") %||% "Covariate model"))
  }
  stop(paste0("wald_test() needs a growth model from fit_gmm() or fit_lcga(), ",
              "or a mixture model with a covariate structural model."),
       call. = FALSE)
}

# The constraint matrix R, one row per constraint, one column per parameter.
.wald_constraint_matrix <- function(constraints, labels) {
  check <- function(nm) {
    bad <- setdiff(nm, labels)
    if (length(bad))
      stop(sprintf("Unknown parameter label%s: %s.\nAvailable labels:\n  %s",
                   if (length(bad) > 1) "s" else "",
                   paste0('"', bad, '"', collapse = ", "),
                   paste(labels, collapse = "\n  ")), call. = FALSE)
  }
  p <- length(labels)
  if (is.matrix(constraints)) {
    if (!is.numeric(constraints) || is.null(colnames(constraints)))
      stop("A constraint matrix must be numeric with parameter labels as column names.",
           call. = FALSE)
    check(colnames(constraints))
    R <- matrix(0, nrow(constraints), p)
    for (j in seq_len(ncol(constraints))) {
      k <- match(colnames(constraints)[j], labels)
      R[, k] <- R[, k] + constraints[, j]
    }
    return(R)
  }
  if (!is.list(constraints) || !length(constraints))
    stop("`constraints` must be a non-empty list or a matrix.", call. = FALSE)
  R <- matrix(0, length(constraints), p)
  for (i in seq_along(constraints)) {
    ci <- constraints[[i]]
    if (is.character(ci)) {
      if (length(ci) != 2L)
        stop("A constraint given as labels must be a pair: c(\"a\", \"b\") for a = b.",
             call. = FALSE)
      ci <- stats::setNames(c(1, -1), ci)
    }
    if (!is.numeric(ci) || is.null(names(ci)) || any(names(ci) == ""))
      stop(paste0("Each constraint must be a pair of labels or a named numeric ",
                  "vector of weights."), call. = FALSE)
    check(names(ci))
    for (j in seq_along(ci)) {
      k <- match(names(ci)[j], labels)
      R[i, k] <- R[i, k] + ci[[j]]
    }
  }
  R
}
