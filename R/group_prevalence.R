# ==============================================================================
# S3 Group-Prevalence Structural Model (a class held equal across groups)
# ==============================================================================
#
# `group_effects = "prevalence"` (or `"both"`) already lets every class
# prevalence vary freely by group, by entering `group` as a class-membership
# covariate through the ordinary multinomial-logit `covariate` model
# (R/covariate.R). That parameterisation has no way to hold *one* class's
# prevalence equal across groups while leaving the others free -- a logit
# coefficient of zero on the group dummies pins that class to have the same
# prevalence as the reference group specifically, not the same prevalence in
# every group.
#
# This emission instead parameterises the per-group class probabilities
# directly, as a `G x K` matrix `gamma`, and lets a caller-chosen subset of
# columns (`frozen`) be held to one shared value across every row while the
# rest renormalise freely within each group. It meets the same contract
# `log_likelihood.covariate()` does -- an `n x K` matrix of
# `log P(class | group_i)` -- so `e_step()` needs no change beyond the
# `.supplies_class_probs()` predicate already covering it (R/em_core.R).
#
# Verified against the reference `em_mglca()` implementation
# (`collins_lanza_ch5.R`, `smoke_pipeline.R`), which checks its own
# unrestricted fit against this package's `covariate` route before trusting
# the restricted ones; see internal/ROADMAP.md Part 17.4 Item B.

# Constructor. `frozen` is a vector of 1-based class indices held equal across
# groups; empty (the default) fits every group's prevalences freely, the
# saturated model that is algebraically identical to the unconstrained
# `covariate` route on group dummies.
group_prevalence_model <- function(n_components, n_groups, frozen = integer(0), ...) {
  state <- list(n_components = n_components, n_groups = n_groups,
               frozen = as.integer(frozen), parameters = list())
  class(state) <- c("group_prevalence", "emission")
  state
}

# Expected counts n_qc: rows are groups, columns are classes. `resp` is n x K
# responsibilities (or initial, unfitted proportions, from init_params());
# `Y` is the group factor's integer codes, one per row of `resp`. Sample
# weights and the `latent` Dirichlet pseudo-cases are folded in here so both
# init_params() and m_step() project onto the same constrained surface from
# the same counts.
#
# The prior mass is split evenly over all K*G cells (`alpha / (K*G)`), the
# same split `.fit_mnl()`'s Dirichlet augmentation applies when its unique
# covariate patterns are exactly the G groups (`w_prior <- alpha / (K * U0)`,
# R/covariate.R) -- which is what makes the unconstrained fit here reproduce
# the `covariate` route's `ll` and `n_parameters()` exactly rather than
# merely approximately.
.group_prevalence_counts <- function(Y, resp, weights, K, G, alpha) {
  w <- if (is.null(weights)) rep(1, nrow(resp)) else weights
  resp_w <- resp * w
  n_qc <- matrix(0, G, K)
  for (q in seq_len(G)) {
    rows <- which(Y == q)
    if (length(rows)) n_qc[q, ] <- colSums(resp_w[rows, , drop = FALSE])
  }
  n_qc + alpha / (K * G)
}

# The constrained maximiser (see the derivation in internal/ROADMAP.md Part
# 17.4 Item B): freezing class c across groups means gamma[, c] is one number
# s_c = A_c / N for every group, while the unfrozen classes stay free within
# each group and renormalise to whatever mass is left. `S = all classes`
# collapses to `s = colSums(n_qc) / N`, i.e. every row equal to the pooled
# share -- the identity the D_ALL benchmark checks.
.group_prevalence_gamma <- function(n_qc, frozen) {
  G <- nrow(n_qc); K <- ncol(n_qc)
  if (!length(frozen)) return(n_qc / rowSums(n_qc))

  N <- sum(n_qc)
  s <- colSums(n_qc[, frozen, drop = FALSE]) / N
  gamma <- matrix(0, G, K)
  gamma[, frozen] <- matrix(rep(s, each = G), nrow = G)
  if (length(frozen) < K) {
    others <- n_qc[, -frozen, drop = FALSE]
    gamma[, -frozen] <- (1 - sum(s)) * others / rowSums(others)
  }
  gamma
}

#' @exportS3Method
init_params.group_prevalence <- function(model_state, X, resp, ...) {
  K <- model_state$n_components
  G <- model_state$n_groups
  # No responsibilities exist yet the first time the main EM loop builds a
  # model state (em_core.R passes `resp = NULL` here, matching
  # init_params.covariate()'s all-zero beta -- the softmax of which is also
  # uniform). Where real initial responsibilities *are* available -- the
  # stepwise 3-step path, and BCH's own init_params() call -- project them
  # onto the constraint surface immediately, the same discipline
  # init_params.blocks() follows for `invariant_items`.
  if (is.null(resp)) {
    model_state$parameters$gamma <- matrix(1 / K, G, K)
    return(model_state)
  }
  alpha <- .bayes_alpha(model_state, "latent")
  n_qc <- .group_prevalence_counts(X, resp, NULL, K, G, alpha)
  model_state$parameters$gamma <- .group_prevalence_gamma(n_qc, model_state$frozen)
  model_state
}

#' @exportS3Method
m_step.group_prevalence <- function(model_state, X, resp, weights = NULL, ...) {
  K <- model_state$n_components
  G <- model_state$n_groups
  alpha <- .bayes_alpha(model_state, "latent")
  n_qc <- .group_prevalence_counts(X, resp, weights, K, G, alpha)
  model_state$parameters$gamma <- .group_prevalence_gamma(n_qc, model_state$frozen)
  model_state
}

#' @exportS3Method
log_likelihood.group_prevalence <- function(model_state, X, ...) {
  log(pmax(model_state$parameters$gamma[X, , drop = FALSE], 1e-300))
}

#' @exportS3Method
n_parameters.group_prevalence <- function(model_state, ...) {
  K <- model_state$n_components
  G <- model_state$n_groups
  nS <- length(model_state$frozen)
  if (nS < K) G * (K - 1) - nS * (G - 1) else K - 1
}

# ------------------------------------------------------------------------------
# Anchoring a constrained fit to an unconstrained one
# ------------------------------------------------------------------------------

# The `G x K` matrix of per-group class probabilities a fitted multiple-group
# model implies, whichever of the two prevalence parameterisations produced it.
# .group_class_sizes() (R/wrapper.R) reads it from here rather than rebuilding
# it inline.
.group_gamma_matrix <- function(object) {
  gi <- object$group_info
  if (is.null(gi) || is.null(object$sm)) return(NULL)
  if (inherits(object$sm, "group_prevalence"))
    return(object$sm$parameters$gamma)
  if (!inherits(object$sm, "covariate")) return(NULL)
  # One dummy-coded row per group level, in the same treatment-contrast coding
  # .lta_group_design() built the fitted regression on (first level =
  # reference), rather than expanding to one row per case only to drop the
  # duplicates -- a case's row depends on nothing but its own group.
  G <- length(gi$levels)
  D <- matrix(0, G, G - 1L)
  if (G > 1L) for (i in 2:G) D[i, i - 1L] <- 1
  softmax_rows(cbind(1, D) %*% t(object$sm$parameters$beta))
}

# `start_from`: run one fit anchored at another fit's solution, with the random
# restarts off. Built only by fit_mixture()'s `group_prevalence_equal` branch.
#
# "Class k has the same prevalence in every group" is a restriction on a *named*
# class, and a mixture's classes are named only by their own parameters. Permute
# the emission's class rows and the columns of `gamma` together and the frozen
# slot lands on a different substantive class at exactly the same likelihood, so
# the restricted model's *global* maximum is one and the same number for every
# k: whichever class is cheapest to freeze. A search reports that number every
# time, and it is not the answer to the question. What the question asks for is
# the maximum in the basin of the unrestricted solution, where class k still
# means what it meant there.
#
# This is therefore the one warm start in the package allowed to *replace* the
# search rather than join it, and the caller has to ask for it by name.
# Everywhere else -- R/group_blocks.R, and Shireman, Steinley and Brusco (2017)
# behind it -- an informed start competes with the random pool on
# log-likelihood, and that competition is the whole safety argument. Here there
# is nothing left to compete for: every basin the random pool could find is the
# same basin relabelled.
#
# Clogg and Goodman (1985, p. 89) single this restriction out as the one needing
# "special treatment", and the reference program reaches the same arrangement
# from the other side: its per-class gamma restriction is written as an
# equivalence set, and supplying starting values (which is how the restriction
# is anchored) makes its random-start option unavailable in the same call.
.group_prevalence_warm_start <- function(donor) {
  force(donor)
  gamma0 <- .group_gamma_matrix(donor)
  function(model_state, X, Y) {
    if (is.null(gamma0) || is.null(donor$mm) || is.null(model_state$sm) ||
        !identical(dim(gamma0),
                   c(as.integer(model_state$sm$n_groups),
                     as.integer(model_state$n_components))))
      stop("`start_from` does not describe the same model as this fit.",
           call. = FALSE)
    # A freshly built emission carries no parameters at all, so it is
    # initialised here for its shape and every value is then overwritten from
    # the donor.
    mm <- .copy_emission_parameters(init_params(model_state$mm, X, NULL),
                                    donor$mm)
    if (is.null(mm))
      stop("`start_from`'s measurement model is not the same shape as this ",
           "fit's; the anchor must differ from it only in ",
           "`group_prevalence_equal`.", call. = FALSE)
    model_state$mm      <- mm
    model_state$weights <- donor$weights
    # Deliberately *not* projected onto the constraint surface, unlike
    # init_params.group_prevalence(). The seed is the unrestricted solution,
    # which is the whole point of anchoring; the first M-step projects it.
    model_state$sm$parameters$gamma <- gamma0
    model_state
  }
}
