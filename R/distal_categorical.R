# ==============================================================================
# S3 Distal Categorical Outcome (No Covariate)
# ==============================================================================
# Estimates class-specific probabilities for a categorical (binary or
# polytomous) distal outcome with no covariate adjustment. This is the
# categorical analogue of distal_continuous: the structural model contains
# only class intercepts, one per outcome category (minus the reference).
#
# All S3 methods (init_params, m_step, log_likelihood, n_parameters) are
# inherited from distal_pooled via S3 class inheritance. distal_pooled handles
# the no-covariate case correctly when Y has exactly one column (D_cov = 0,
# L = K), so the design matrix reduces to a K-column class-indicator matrix
# and beta_pooled holds only class intercepts.
#
# Recommended correction: ML (same as distal_pooled / distal_regression).
# ==============================================================================

distal_categorical_model <- function(n_components, tol = 1e-4, max_iter = 500,
                                     method = "newton-raphson") {
  state <- list(n_components = n_components, tol = tol, max_iter = max_iter,
                method = method, parameters = list())
  class(state) <- c("distal_categorical", "distal_pooled", "emission")
  return(state)
}
