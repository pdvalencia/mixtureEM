# ==============================================================================
# Class-vs-class contrasts on a distal outcome
# ==============================================================================
# The omnibus test printed by summary() answers "do the classes differ on this
# outcome at all", which is the question a reviewer asks first and almost never
# the question the paper is about. What gets written up is which classes differ:
# the high-risk class scores half a standard deviation above the normative one,
# the two intermediate classes are indistinguishable. Those are contrasts, and
# every one of them is a linear combination of parameters the fit already
# carries, together with a covariance matrix it already estimated.
#
# The covariance is the whole reason this belongs in the package rather than in
# a user's script. Two class means from the same fit are not independent - the
# same posteriors, and under BCH the same inverted classification-error matrix,
# enter both - so the standard error of their difference is not the root of the
# sum of their squared standard errors. Off the diagonal of Sigma_mu the two
# can differ by a good deal, and always in the direction that matters: positive
# covariance between the two means makes the naive interval too wide and the
# test too conservative, so the difference a user computes by hand is the one
# they are most likely to miss.

# Find a distal sub-model of a given class, whether the structural model is that
# model or a nested one holding it alongside a covariate model. Three copies of
# this unwrapping already existed in summary(); this is the one the new code
# uses, and `exclude` covers the one case where the classes overlap
# (distal_categorical inherits from distal_pooled).
.distal_submodel <- function(sm, what, exclude = NULL) {
  ok <- function(m) {
    if (is.null(m) || !inherits(m, what)) return(FALSE)
    if (!is.null(exclude) && inherits(m, exclude)) return(FALSE)
    TRUE
  }
  if (ok(sm)) return(sm)
  if (inherits(sm, "nested") && "distal" %in% names(sm$models) &&
      ok(sm$models$distal))
    return(sm$models$distal)
  NULL
}

# The distal blocks of a nested structural model, named by the outcome each one
# models (the first column of its block of `Y`). A joint model of predictors
# and several outcomes holds them as "distal", "distal2", ...
.distal_blocks <- function(sm, Y) {
  if (!inherits(sm, "nested")) return(character(0))
  ends <- cumsum(sm$columns_per_model)
  nm   <- names(sm$models)
  keep <- startsWith(nm, "distal")
  first <- (ends - sm$columns_per_model + 1)[keep]
  lab  <- colnames(Y)[first] %||% nm[keep]
  stats::setNames(nm[keep], lab)
}

# A view of `fit` in which the requested outcome's block is the one named
# "distal", so every reader written for a single distal outcome serves any of
# several. `outcome` is an outcome label or a block position; NULL is the first.
.select_distal <- function(fit, outcome = NULL) {
  blocks <- .distal_blocks(fit$sm, fit$Y)
  if (length(blocks) < 2L) {
    if (!is.null(outcome) && length(blocks) == 1L &&
        !outcome %in% c(names(blocks), 1))
      stop(sprintf("`outcome` \"%s\" is not on this fit; its outcome is %s.",
                   outcome, names(blocks)), call. = FALSE)
    return(fit)
  }
  pick <- if (is.null(outcome)) 1L
          else if (is.numeric(outcome)) as.integer(outcome)
          else match(outcome, names(blocks))
  if (length(pick) != 1L || is.na(pick) || pick < 1L || pick > length(blocks))
    stop(sprintf("`outcome` must name one of this fit's outcomes: %s.",
                 paste(names(blocks), collapse = ", ")), call. = FALSE)
  fit$sm$models$distal <- fit$sm$models[[blocks[[pick]]]]
  attr(fit, "distal_label") <- names(blocks)[pick]
  fit
}

# The (class, reference) pairs to report, in printing order.
#
# With no reference given, every unordered pair, each reported once as the
# higher-numbered class against the lower. With one given, that class is held
# fixed on the right-hand side, which is the layout the categorical outcome
# tables already use.
.contrast_pairs <- function(K, ref = NULL) {
  if (is.null(ref)) {
    pairs <- utils::combn(K, 2L)
    return(data.frame(class = pairs[2L, ], reference = pairs[1L, ]))
  }
  data.frame(class = setdiff(seq_len(K), ref), reference = ref)
}

# One contrast: the difference between two entries of theta, and its standard
# error from the corresponding 2x2 block of V.
#
# Var(theta_a - theta_b) = V[a,a] + V[b,b] - 2 V[a,b]. The cross term is the
# point; dropping it is the mistake this function exists to prevent.
.contrast_row <- function(theta, V, a, b) {
  est <- theta[a] - theta[b]
  var <- V[a, a] + V[b, b] - 2 * V[a, b]
  se  <- sqrt(max(0, var))
  z   <- if (is.finite(se) && se > 0) est / se else NA_real_
  p   <- if (is.na(z)) NA_real_ else 2 * stats::pnorm(-abs(z))
  list(estimate = est, se = se, z = z, p = p)
}

#' Class-vs-Class Contrasts on a Distal Outcome
#'
#' @description
#' Tests which latent classes differ on a distal outcome attached by
#' [add_outcome()], rather than whether any of them do. The omnibus Wald test
#' printed by [summary()] is a single joint statistic; this is the set of
#' contrasts behind it, each with a standard error, a confidence interval and a
#' p-value.
#'
#' Each contrast is reported as `class` minus `reference`, so a positive
#' estimate means the class scores above the class it is compared with. With
#' `ref = NULL` every pair of classes is reported once; naming a class in `ref`
#' holds it fixed as the comparison.
#'
#' The standard errors account for the covariance between the class parameters,
#' which is not optional here. Two class means from one fit are estimated from
#' the same posteriors, and under the BCH correction from the same inverted
#' classification-error matrix, so they are correlated; the standard error of
#' their difference is not the root of the sum of their squared standard errors,
#' and the naive version is conservative exactly when the correlation is
#' positive. The full sandwich covariance the fit already carries is used
#' wherever it is available.
#'
#' @param fit A fitted model with a distal outcome attached, as returned by
#'   [add_outcome()], or a three-step [fit_lta()] fit with `distal`, where the
#'   classes are the statuses at the last occasion.
#' @param ref Optional reference class, as an integer between 1 and the number
#'   of classes. `NULL` (the default) reports all pairs.
#' @param adjust Multiplicity adjustment across the reported contrasts, passed
#'   to [stats::p.adjust()]: `"none"` (the default), `"BH"`, `"holm"` or
#'   `"bonferroni"`. With `"none"` the p-values are the per-contrast ones, which
#'   is what to report when the contrasts were named in advance; an adjustment
#'   belongs on the all-pairs table read as a family. See "Choosing an
#'   adjustment" below.
#' @param level Confidence level for the intervals. Default `0.95`.
#' @param outcome For a fit holding several distal outcomes (see `predictors`
#'   in [add_outcome()]), the outcome to contrast, by name or position.
#'   `NULL` (the default) is the first; the `outcome` attribute of the result
#'   names the one used.
#' @param term For an outcome regressed on covariates with class-specific
#'   slopes (`slopes` in [add_outcome()] naming the covariate), the covariate
#'   whose slopes to contrast, by its column name. `NULL` (the default)
#'   contrasts the class intercepts instead. See "Slope contrasts" below.
#' @param effect_size If `TRUE`, add standardized effect sizes for a
#'   continuous outcome: Cohen's d for every contrast and the LTB-\eqn{\omega}
#'   across all classes. Default `FALSE`. See "Effect sizes" below.
#' @param ... Currently unused.
#'
#' @section Effect sizes:
#' A contrast's p-value says whether two classes can be told apart, not how
#' far apart they are; with enough cases a trivial gap is significant, and a
#' large one can miss significance in a small class. `effect_size = TRUE`
#' reports the size of the differences on a scale that does not depend on the
#' outcome's units.
#'
#' For each pair, Cohen's d is the difference between the two class means
#' divided by the square root of the average of the two class variances,
#' \eqn{d_{jk} = (\mu_j - \mu_k) / \sqrt{(\sigma^2_j + \sigma^2_k)/2}}
#' (Carter, 2026); with a variance shared by the classes it is the gap
#' measured in within-class standard deviations. The point estimate is the
#' standardized mean difference of Hedges and Olkin (2014).
#'
#' Across all classes, the LTB-\eqn{\omega} of Lanza, Tan and Bray (2013) is
#' the spread of the class means around their overall mean, each weighted by
#' its class size: \eqn{\omega = \sqrt{\sum_k \pi_k (\mu_k - \bar\mu)^2}},
#' \eqn{\bar\mu = \sum_k \pi_k \mu_k}. It is in the outcome's units. The
#' printed value is standardized by the outcome's total standard deviation,
#' \eqn{\omega_{std} = \omega / \sqrt{\sum_k \pi_k \sigma^2_k + \omega^2}}, so
#' that \eqn{\omega_{std}^2} is the share of the outcome's variance that lies
#' between the classes. Both are in the `omega` attribute of the result
#' (`$standardized` and `$raw`, each with `estimate`, `se`, `lower` and
#' `upper`); the raw one is the \eqn{\omega} of Carter's reporting workflow.
#' The effect size answers a different question from the omnibus test in
#' [summary()] and does not replace it.
#'
#' Standard errors and intervals are by the delta method from the joint
#' covariance of the class sizes, means and variances, so the uncertainty in
#' the sizes and the variances is carried along with that of the means. That
#' covariance is stored by a one-step [fit_gmm()] or [fit_lcga()] outcome
#' (the robust sandwich of [growth_estimates()]) and by
#' `add_outcome(correction = "ML")`; for other outcome models, including the
#' BCH correction, `effect_size = TRUE` is an error. No thresholds such as
#' "small" or "large" are attached: what size of difference matters depends on
#' the outcome and the field.
#'
#' @section Choosing an adjustment:
#' The adjustments answer different questions. `"holm"` and `"bonferroni"`
#' control the family-wise error rate, the chance of even one false positive
#' among the contrasts; that is the right guarantee for a narrowly
#' confirmatory comparison, where a single wrong claim is costly. `"BH"`
#' (Benjamini and Hochberg, 1995) controls the false discovery rate, the
#' expected share of false positives among the contrasts declared
#' significant; that is the right guarantee when the aim is to describe the
#' pattern of differences between the classes, which is the usual aim of a
#' distal-outcome analysis.
#'
#' The pairwise contrasts are positively dependent by construction: each class
#' mean enters several of them (in a three-class model, class 1's mean is in
#' both the 1-vs-2 and the 1-vs-3 contrast). Family-wise procedures do not
#' use that dependence and are conservative here, so they suppress more of a
#' real pattern than the number of contrasts alone suggests. Benjamini and
#' Yekutieli (2001) prove that the BH procedure keeps its false discovery rate
#' under positive regression dependence, which covers positively correlated
#' one-sided tests; the p-values here are two-sided, for which that proof does
#' not apply as it stands.
#'
#' @section Slope contrasts:
#' When the class moderates a covariate's effect on the outcome -- each class
#' has its own regression line -- the question after "do the slopes differ at
#' all" (the Wald test printed by [summary()]) is "which classes differ in
#' slope". A slope contrast answers it: the estimate is how much more the
#' outcome rises per unit of the covariate in `class` than in `reference`.
#' This is the second stage of the test sequence for latent class moderation
#' in Arch, Nylund-Gibson and Ing (2026); the first is the omnibus test and the
#' third is each class's slope against zero, both already in [summary()].
#'
#' Fit the outcome with `slopes` naming the covariate, e.g.
#' `add_outcome(fit, y, covariates = df["x"], slopes = "x")`, and pass
#' `term = "x"`. That fit estimates every class's line in one regression, so
#' its covariance includes the covariance between the slopes of different
#' classes, which the contrast standard error needs. On a continuous outcome
#' `slopes = "class_specific"` is the same fit, with every covariate
#' moderated, and can be contrasted the same way.
#'
#' With `term` left `NULL` on such a fit the contrasts are between the class
#' intercepts, i.e. the class means at covariate values of zero; centring the
#' covariate first moves that comparison to its mean.
#'
#' @references
#' Arch, D. A. N., Nylund-Gibson, K., & Ing, M. (2026). Moderation with a
#' latent class variable. \emph{Behavior Research Methods}, 58, 108.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#' A practical and powerful approach to multiple testing. \emph{Journal of the
#' Royal Statistical Society: Series B}, \emph{57}(1), 289–300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery
#' rate in multiple testing under dependency. \emph{The Annals of Statistics},
#' \emph{29}(4), 1165–1188. \doi{10.1214/aos/1013699998}
#'
#' Carter, D. B. (2026). Distal outcomes in mixture modeling: A guide for
#' pairwise comparisons, multiplicity control, and effect size reporting.
#' \emph{Behavior Research Methods}, \emph{58}, 268.
#' \doi{10.3758/s13428-026-03144-4}
#'
#' Hedges, L. V., & Olkin, I. (2014). \emph{Statistical methods for
#' meta-analysis}. Academic Press.
#'
#' Lanza, S. T., Tan, X., & Bray, B. C. (2013). Latent class analysis with
#' distal outcomes: A flexible model-based approach. \emph{Structural Equation
#' Modeling}, \emph{20}(1), 1–26. \doi{10.1080/10705511.2013.742377}
#'
#' @return A data frame of class `outcome_contrasts`, one row per contrast, with
#'   columns `category` (the outcome category the contrast is on, `NA` for a
#'   continuous outcome), `class`, `reference`, `estimate`, `se`, `lower`,
#'   `upper`, `z` and `p`, plus `p_adj` when `adjust` is not `"none"` and
#'   `OR`, `OR_lower`, `OR_upper` for a categorical outcome, whose estimates are
#'   log odds. With `effect_size = TRUE`, also `d`, `d_se`, `d_lower` and
#'   `d_upper`, and the `omega` attribute. The `method` attribute names the
#'   covariance the standard errors came from.
#'
#' @seealso [add_outcome()] for fitting the outcome model and [summary()] for
#'   the omnibus test; [wald_omnibus_test()] and [bootstrap_covariates()] for
#'   the bootstrap route, which is what to use for an outcome whose parameters
#'   are estimated one class at a time.
#'
#' @examples
#' set.seed(1)
#' items <- matrix(rbinom(300, 1, 0.5), nrow = 100)
#' bmi   <- rnorm(100, mean = 25)
#' fit   <- fit_mixture(items, n_classes = 2, measurement = "binary")
#' fit_out <- add_outcome(fit, bmi)
#' outcome_contrasts(fit_out)
#'
#' # Which classes differ in the slope of the outcome on a covariate
#' x     <- rnorm(100)
#' fit_m <- add_outcome(fit, bmi + x, covariates = data.frame(x = x),
#'                      slopes = "x")
#' outcome_contrasts(fit_m, term = "x")
#'
#' @export
outcome_contrasts <- function(fit, ref = NULL,
                              adjust = c("none", "BH", "holm", "bonferroni"),
                              level = 0.95, outcome = NULL, term = NULL,
                              effect_size = FALSE, ...) {
  adjust <- match.arg(adjust)
  if (!is.logical(effect_size) || length(effect_size) != 1L ||
      is.na(effect_size))
    stop("`effect_size` must be TRUE or FALSE.", call. = FALSE)
  if (!is.null(term) && (!is.character(term) || length(term) != 1L))
    stop("`term` must be a single covariate name.", call. = FALSE)
  lta <- inherits(fit, "lta_model")
  if (!inherits(fit, "mixture_model") && !lta)
    stop("`fit` must be a fitted mixture model.", call. = FALSE)
  if (lta) {
    if (!is.null(term))
      stop("`term` applies to an outcome with class-specific slopes; a ",
           "fit_lta() distal outcome has none.", call. = FALSE)
    built <- .lta_distal_contrast_parts(fit, outcome)
    attr(fit, "distal_label") <- built$label
    K <- fit$n_statuses
  } else {
    fit <- .select_distal(fit, outcome)
    K <- fit$n_components
  }
  if (is.null(K) || K < 2L)
    stop("Class contrasts need at least two classes.", call. = FALSE)
  if (!is.null(ref)) {
    if (!is.numeric(ref) || length(ref) != 1L || ref < 1 || ref > K)
      stop(sprintf("`ref` must be a class index between 1 and %d.", K),
           call. = FALSE)
    ref <- as.integer(ref)
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1)
    stop("`level` must be a confidence level strictly between 0 and 1.",
         call. = FALSE)

  if (!lta) built <- .outcome_contrast_parts(fit, K, term)
  pairs <- .contrast_pairs(K, ref)

  rows <- list()
  for (part in built$parts) {
    for (i in seq_len(nrow(pairs))) {
      a <- pairs$class[i]; b <- pairs$reference[i]
      r <- .contrast_row(part$theta, part$V, part$index[a], part$index[b])
      rows[[length(rows) + 1L]] <- data.frame(
        category  = part$category,
        class     = a,
        reference = b,
        estimate  = r$estimate, se = r$se, z = r$z, p = r$p,
        stringsAsFactors = FALSE)
    }
  }

  out <- .rbind_tidy(rows)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out$lower <- out$estimate - crit * out$se
  out$upper <- out$estimate + crit * out$se

  if (adjust != "none") out$p_adj <- stats::p.adjust(out$p, method = adjust)

  cols <- c("category", "class", "reference", "estimate", "se",
            "lower", "upper", "z", "p")
  if (adjust != "none") cols <- c(cols, "p_adj")
  out <- out[, cols, drop = FALSE]

  omega <- NULL
  if (effect_size) {
    if (lta || !identical(built$type, "continuous") || !is.null(term))
      stop("Effect sizes are defined for the class means of a continuous ",
           "outcome; this contrast is on ",
           if (!is.null(term)) "a slope." else "another kind of parameter.",
           call. = FALSE)
    inp <- .effect_size_inputs(fit, K)
    es  <- .distal_effect_sizes(inp$theta, inp$V, K, pairs, level)
    out$d <- es$d; out$d_se <- es$d_se
    out$d_lower <- es$d_lower; out$d_upper <- es$d_upper
    omega <- es$omega
  }

  if (identical(built$type, "categorical")) {
    out$OR       <- exp(out$estimate)
    out$OR_lower <- exp(out$lower)
    out$OR_upper <- exp(out$upper)
  }

  rownames(out) <- NULL
  attr(out, "method")       <- built$method
  attr(out, "outcome_type") <- built$type
  attr(out, "level")        <- level
  attr(out, "adjust")       <- adjust
  attr(out, "outcome")      <- attr(fit, "distal_label")
  attr(out, "omega")        <- omega
  class(out) <- c("outcome_contrasts", "data.frame")
  out
}

# Reduce whichever distal model is attached to the pieces a contrast needs: a
# parameter vector, its covariance, and the positions in that covariance the
# classes occupy. One "part" per outcome category, so a polytomous outcome is a
# stack of the same contrasts on each of its non-reference categories.
.outcome_contrast_parts <- function(fit, K, term = NULL) {

  cpool <- .distal_submodel(fit$sm, "distal_continuous_pooled")
  if (!is.null(term)) {
    if (is.null(cpool) || is.null(cpool$parameters$cov_theta))
      stop("`term` names a class-specific slope, and this outcome model has ",
           "none. Fit the outcome with `slopes` naming the covariate.",
           call. = FALSE)
    return(.slope_contrast_parts(cpool, K, term))
  }

  # --- Continuous outcome: class means -----------------------------------
  cont <- .distal_submodel(fit$sm, "distal_continuous")
  if (!is.null(cont) && !is.null(cont$parameters$means)) {
    means <- as.vector(cont$parameters$means)
    Sigma <- cont$parameters$Sigma_mu

    # Sigma_mu is stored only by the BCH correction. Without it the class means
    # are treated as independent, which is what the omnibus test falls back to
    # as well; the attribute says so, because a user comparing the two needs to
    # know which one they got.
    if (is.null(Sigma) || !all(is.finite(Sigma))) {
      ses    <- as.vector(cont$parameters$ses)
      Sigma  <- diag(ses^2, nrow = K)
      method <- "independent class means (no cross-class covariance stored)"
    } else {
      method <- "sandwich covariance of the class means"
    }

    return(list(
      type   = "continuous",
      method = method,
      parts  = list(list(category = NA_integer_, theta = means, V = Sigma,
                         index = seq_len(K)))))
  }

  # --- Continuous outcome adjusted for covariates: class intercepts --------
  #
  # The first K entries of beta_pooled are the class intercepts, and cov_theta
  # is their covariance with the slopes; the difference of two intercepts is
  # the difference of the class means at equal covariate values.
  if (!is.null(cpool) && !is.null(cpool$parameters$cov_theta)) {
    return(list(
      type   = "continuous",
      method = "covariance of the class intercepts (class means at covariates of zero)",
      parts  = list(list(category = NA_integer_,
                         theta = as.vector(cpool$parameters$beta_pooled)[seq_len(K)],
                         V = cpool$parameters$cov_theta[seq_len(K), seq_len(K),
                                                        drop = FALSE],
                         index = seq_len(K)))))
  }

  # --- Categorical outcome: class intercepts, one set per category --------
  #
  # beta_pooled is (M-1) x L with L = K + D_cov, and is vectorised row-major
  # into the Hessian, so category m's intercept for class k sits at
  # (m-1)*L + k. That indexing is fixed by .wald_omnibus_pooled() and by the
  # odds-ratio tables in summary(); it is repeated rather than re-derived.
  pooled <- .distal_submodel(fit$sm, "distal_pooled")
  if (!is.null(pooled) && !is.null(pooled$parameters$beta_pooled) &&
      nrow(pooled$parameters$beta_pooled) > 0L) {
    beta <- pooled$parameters$beta_pooled
    if (is.null(pooled$parameters$hessian))
      stop("This outcome model carries no covariance for its parameters, so ",
           "contrasts cannot be formed. Refit it with `se = \"corrected\"` or ",
           "`se = \"robust\"`.", call. = FALSE)

    V <- pinv(-pooled$parameters$hessian)
    L <- ncol(beta)

    parts <- lapply(seq_len(nrow(beta)), function(m) {
      list(category = m + 1L,
           theta    = as.vector(t(beta))[((m - 1L) * L + 1L):(m * L)],
           V        = V[((m - 1L) * L + 1L):(m * L),
                        ((m - 1L) * L + 1L):(m * L), drop = FALSE],
           index    = seq_len(K))
    })

    return(list(type = "categorical",
                method = "robust sandwich covariance of the outcome model",
                parts = parts))
  }

  # --- Everything else ----------------------------------------------------
  # slopes = "class_specific" fits each class's regression on its own and
  # keeps no covariance between classes. The same model fitted with every
  # covariate named in `slopes` is one stacked regression whose sandwich
  # carries that covariance, so that is where a continuous outcome is sent.
  if (!is.null(.distal_submodel(fit$sm, "distal_continuous_regression")))
    stop("This outcome model estimates each class's regression separately ",
         "and stores no covariance between the classes, so class-vs-class ",
         "contrasts cannot be formed from it. Refit with `slopes` naming ",
         "every covariate (e.g. `slopes = c(\"x1\", \"x2\")`), which gives ",
         "the same estimates and keeps that covariance; then use `term` ",
         "to contrast a slope.", call. = FALSE)
  if (!is.null(.distal_submodel(fit$sm, "distal_regression")))
    stop("This outcome model estimates a separate parameter block per class ",
         "and stores no covariance between the blocks, so class-vs-class ",
         "contrasts cannot be formed from it. Fit the outcome with ",
         "`slopes = \"pooled\"`, or use bootstrap_covariates() with ",
         "wald_omnibus_test().", call. = FALSE)

  stop("No distal outcome found on this fit. Attach one with add_outcome() ",
       "first.", call. = FALSE)
}

# The class-specific slopes on one covariate: a block of K entries of
# beta_pooled, located by the "term:Class{k}" column names R/stepwise.R gives
# it, and the matching K x K block of cov_theta. The block is taken from the
# joint covariance, so the covariance between two classes' slopes -- which
# the stacked regression estimates -- enters every contrast.
.slope_contrast_parts <- function(cpool, K, term) {
  beta <- cpool$parameters$beta_pooled
  nm   <- colnames(beta)
  idx  <- match(paste0(term, ":Class", seq_len(K)), nm)
  if (anyNA(idx)) {
    mods <- unique(sub(":Class[0-9]+$", "",
                       grep(":Class[0-9]+$", nm, value = TRUE)))
    stop(sprintf(paste0("`term` \"%s\" has no class-specific slopes on this ",
                        "fit. Class-specific: %s. A covariate with one slope ",
                        "shared by every class has nothing to contrast."),
                 term, if (length(mods)) paste(mods, collapse = ", ")
                       else "(none)"),
         call. = FALSE)
  }
  list(
    type   = "continuous",
    method = sprintf("covariance of the class-specific slopes on %s", term),
    parts  = list(list(category = NA_integer_,
                       theta = as.vector(beta)[idx],
                       V     = cpool$parameters$cov_theta[idx, idx, drop = FALSE],
                       index = seq_len(K))))
}

# Effect sizes on a continuous outcome, from theta = (pi_1..pi_K, mu_1..mu_K,
# V_1..V_K) and its 3K x 3K covariance. That layout, rather than the size
# logits, is what the fits store, because re-sorting the classes then permutes
# each block of K and nothing else. A variance shared by the classes appears K
# times with identical rows, which the delta method below needs no help with.
#
# Pairwise d divides by the average of the two class variances; omega is the
# pi-weighted spread of the class means around their weighted mean, in outcome
# units, and omega_std is omega over the total outcome SD, so omega_std^2 is the
# share of the outcome's variance that lies between the classes.
.effect_size_values <- function(theta, K, pairs) {
  p  <- theta[seq_len(K)]
  mu <- theta[K + seq_len(K)]
  v  <- theta[2L * K + seq_len(K)]
  a  <- pairs$class; b <- pairs$reference
  d  <- (mu[a] - mu[b]) / sqrt((v[a] + v[b]) / 2)
  om <- sqrt(sum(p * (mu - sum(p * mu))^2))
  c(d, omega = om, omega_std = om / sqrt(sum(p * v) + om^2))
}

.distal_effect_sizes <- function(theta, V, K, pairs, level) {
  f   <- function(x) .effect_size_values(x, K, pairs)
  est <- f(theta)
  h   <- 1e-6 * pmax(1, abs(theta))
  G   <- vapply(seq_along(theta), function(j) {
    a <- theta; a[j] <- a[j] + h[j]
    b <- theta; b[j] <- b[j] - h[j]
    (f(a) - f(b)) / (2 * h[j])
  }, numeric(length(est)))
  G   <- matrix(G, length(est))
  se  <- sqrt(pmax(rowSums((G %*% V) * G), 0))
  crit <- stats::qnorm(1 - (1 - level) / 2)
  one <- function(i) list(estimate = est[[i]], se = se[i],
                          lower = est[[i]] - crit * se[i],
                          upper = est[[i]] + crit * se[i])
  np  <- nrow(pairs)
  list(d = est[seq_len(np)], d_se = se[seq_len(np)],
       d_lower = est[seq_len(np)] - crit * se[seq_len(np)],
       d_upper = est[seq_len(np)] + crit * se[seq_len(np)],
       omega = list(standardized = one(np + 2L), raw = one(np + 1L)))
}

# theta and its covariance for the fits that carry one: a one-step growth
# model's joint sandwich (growth_estimates(se = "robust")), carried from the
# size logits to the sizes, or the ML step-3 matrix stored as V_pmv. Anything
# else stops, saying which fits have it.
.effect_size_inputs <- function(fit, K) {
  if (inherits(fit, c("gmm", "lcga")) &&
      identical(as.integer(fit$n_steps), 1L)) {
    ge <- growth_estimates(fit, se = "robust")
    Vg <- attr(ge, "vcov")
    il <- which(ge$parameter == "class_logit")
    im <- which(ge$parameter == "outcome_mean")
    iv <- which(ge$parameter == "outcome_variance")
    if (!is.null(Vg) && length(il) == K - 1L && length(im) == K &&
        length(iv) == K) {
      a  <- c(ge$estimate[il], 0)
      p  <- exp(a - max(a)); p <- p / sum(p)
      # d pi / d logit_j = pi * (1[k = j] - pi_j); the class-K logit is fixed.
      Jp <- (diag(p) - outer(p, p))[, seq_len(K - 1L), drop = FALSE]
      idx <- c(il, im, iv)
      J <- matrix(0, 3L * K, length(idx))
      J[seq_len(K), seq_len(K - 1L)] <- Jp
      J[K + seq_len(2L * K), K - 1L + seq_len(2L * K)] <- diag(2L * K)
      return(list(theta = c(p, ge$estimate[im], ge$estimate[iv]),
                  V = J %*% Vg[idx, idx, drop = FALSE] %*% t(J)))
    }
  }
  cont <- .distal_submodel(fit$sm, "distal_continuous")
  if (!is.null(cont) && !is.null(cont$parameters$V_pmv)) {
    p <- cont$parameters
    return(list(theta = c(p$pk_step3, p$means[, 1L],
                          rep_len(p$covariances[, 1L], K)),
                V = p$V_pmv))
  }
  stop("Effect sizes need the joint covariance of the class sizes, means and ",
       "variances of a continuous outcome, which only a one-step fit_gmm() or ",
       "fit_lcga() outcome and an add_outcome(correction = \"ML\") outcome ",
       "carry. They are not available for this outcome model.", call. = FALSE)
}

#' @export
print.outcome_contrasts <- function(x, ...) {
  cat("Class contrasts on the distal outcome\n")
  cat(sprintf("Standard errors: %s\n", attr(x, "method")))
  if (!identical(attr(x, "adjust"), "none"))
    cat(sprintf("P-values adjusted across %d contrasts (%s)\n",
                nrow(x), attr(x, "adjust")))
  om <- attr(x, "omega")$standardized
  if (!is.null(om))
    cat(sprintf("Standardized LTB-omega = %.3f [%.3f, %.3f]\n",
                om$estimate, om$lower, om$upper))
  cat("\n")

  poly <- any(!is.na(x$category))
  lvl  <- sprintf("[%g%% CI]", 100 * attr(x, "level"))
  # The widths here track the value format below exactly; the interval field is
  # 18 characters wide ("[" + 7 + ", " + 7 + "]"), so its heading is padded to
  # the same and the columns line up.
  cat(sprintf("  %-13s %9s %9s  %-18s  %s\n",
              if (poly) "Cat  Contrast" else "Contrast",
              "Estimate", "SE", lvl, "P-Value"))

  p_col <- if ("p_adj" %in% names(x)) x$p_adj else x$p
  for (i in seq_len(nrow(x))) {
    lab <- sprintf("%d vs %d", x$class[i], x$reference[i])
    if (poly) lab <- sprintf("%-4d %s", x$category[i], lab)
    cat(sprintf("  %-13s %9.3f %9.3f  [%7.3f, %7.3f]  %s\n",
                lab, x$estimate[i], x$se[i], x$lower[i], x$upper[i],
                .fmt_pval(p_col[i])))
  }
  if ("d" %in% names(x)) {
    cat(sprintf("\n  %-13s %9s %9s  %-18s\n", "Contrast", "d", "SE", lvl))
    for (i in seq_len(nrow(x)))
      cat(sprintf("  %-13s %9.3f %9.3f  [%7.3f, %7.3f]\n",
                  sprintf("%d vs %d", x$class[i], x$reference[i]),
                  x$d[i], x$d_se[i], x$d_lower[i], x$d_upper[i]))
  }
  invisible(x)
}
