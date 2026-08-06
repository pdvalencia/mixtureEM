# ==============================================================================
# Bayes constants (prior strengths)
# ==============================================================================
#
# Every M-step in the package maximises a *penalised* likelihood: a weak prior
# is added to each block of parameters so that a class supported by very little
# data is pulled towards a conservative null rather than towards a boundary. The
# strengths of those priors used to be hard-coded at 1. They are now one list,
# carried on the model state, so that a user reproducing an unregularized
# reference fit can switch them off and a user whose fit has collapsed can turn
# the relevant one up.
#
# The four names, and what each one guards:
#
#   latent      Dirichlet prior on the class weights (and, in the transition
#               models, on the initial-status and transition probabilities).
#               Stops a class from being emptied outright.
#   categorical Dirichlet prior on item-response probabilities, centred on the
#               item's observed marginal. Stops a response probability from
#               being pinned at 0 or 1 by a handful of cases.
#   poisson     The Gamma analogue of the above for count indicators.
#   variances   Truncated inverse-Wishart prior on the class-specific variances
#               of continuous indicators, centred on the item's observed
#               marginal variance. Stops a variance from collapsing to zero.
#
# Each is expressed in *pseudo-observations*: a constant of alpha contributes
# alpha / K pseudo-observations to each of the K classes, so its influence
# shrinks as the sample grows and it is negligible on any class with real data
# behind it. Setting one to 0 removes that prior entirely and recovers plain
# maximum likelihood for that block.
#
# On why the variance prior is written in the *truncated* inverse-Wishart form
# (no -((K+1)/2) log|Sigma| term), see the comment in m_step.gaussian_diag().
# The short version: the truncated form is the one for which alpha = 0 recovers
# ML exactly, which is the whole point of the escape hatch.

# The defaults. All four at 1, which is the value the first three were already
# hard-coded to; `variances = 1` is new and is deliberately weak — see the
# `bayes_constants` documentation in fit_mixture() for why the stronger value
# that fixes a collapsed fit is a recommendation rather than a default.
.bayes_constants_defaults <- list(
  latent      = 1,
  categorical = 1,
  poisson     = 1,
  variances   = 1
)

# Validate a user-supplied `bayes_constants` and fill in the defaults.
.resolve_bayes_constants <- function(bc = NULL) {
  if (is.null(bc)) return(.bayes_constants_defaults)
  if (!is.list(bc))
    stop("`bayes_constants` must be a list, e.g. list(variances = 5).",
         call. = FALSE)

  unknown <- setdiff(names(bc), names(.bayes_constants_defaults))
  if (length(unknown))
    stop(sprintf("Unknown `bayes_constants`: %s. Valid names are %s.",
                 paste(sQuote(unknown), collapse = ", "),
                 paste(names(.bayes_constants_defaults), collapse = ", ")),
         call. = FALSE)

  out <- .bayes_constants_defaults
  for (nm in names(bc)) {
    v <- bc[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0)
      stop(sprintf("`bayes_constants$%s` must be a single non-negative number.",
                   nm), call. = FALSE)
    out[[nm]] <- as.numeric(v)
  }
  out
}

# Attach the resolved constants to an emission and, recursively, to every
# sub-model inside it.
#
# Recursion is what makes "the M-step and the polish read the same object"
# literally true for the block and nested models: `m_step.blocks()` dispatches
# to a sub-model, and `.refine_time_block_view()` builds its flat view out of
# the first sub-model, so both arrive at the same attached list. An emission
# that carries no continuous parameters simply ignores it.
.attach_bayes_constants <- function(mm, bc) {
  if (is.null(mm)) return(NULL)
  mm$bayes_constants <- bc
  if (!is.null(mm$models) && is.list(mm$models))
    mm$models <- lapply(mm$models, .attach_bayes_constants, bc = bc)
  mm
}

# Read one constant off an emission (or any object carrying the list), falling
# back to the package default when the object predates the plumbing — a fitted
# object reloaded from an older session, or an emission constructed directly in
# a test rather than through fit_mixture().
.bayes_alpha <- function(x, name) {
  v <- x$bayes_constants[[name]]
  if (is.null(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0)
    return(.bayes_constants_defaults[[name]])
  v
}
