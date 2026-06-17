# ==============================================================================
# S3 Emission Factory (List Dispatcher)
# ==============================================================================

# Build an emission model state based on a string or list descriptor.
build_emission <- function(descriptor, n_components = 2, ...) {

  # 1. Handle Nested Models (List of sub-models)
  if (is.list(descriptor)) {
    state <- list(
      n_components      = n_components,
      models            = list(),
      columns_per_model = numeric()
    )
    class(state) <- c("nested", "emission")

    for (name in names(descriptor)) {
      item       <- descriptor[[name]]
      model_type <- item$model
      n_cols     <- item$n_columns

      args            <- item
      args$model      <- NULL
      args$n_columns  <- NULL
      args$n_components <- n_components

      state$models[[name]] <- do.call(build_emission,
                                      c(list(descriptor = model_type), args))
      state$columns_per_model[name] <- n_cols
    }
    return(state)
  }

  # 2. Handle Simple Models (Character string)
  if (is.character(descriptor)) {
    args_list <- list(n_components = n_components, ...)

    # ------------------------------------------------------------------
    # Measurement models
    # ------------------------------------------------------------------
    if (descriptor %in% c("bernoulli", "binary")) {
      return(do.call(categorical_model, c(list(type = "bernoulli"), args_list)))
    } else if (descriptor %in% c("bernoulli_nan", "binary_nan")) {
      return(do.call(categorical_model, c(list(type = "bernoulli_nan"), args_list)))
    } else if (descriptor %in% c("multinoulli", "categorical")) {
      return(do.call(categorical_model, c(list(type = "multinoulli"), args_list)))
    } else if (descriptor %in% c("multinoulli_nan", "categorical_nan")) {
      return(do.call(categorical_model, c(list(type = "multinoulli_nan"), args_list)))
    } else if (descriptor %in% c("gaussian_unit", "gaussian")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_unit"), args_list)))
    } else if (descriptor %in% c("gaussian_unit_nan", "gaussian_nan")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_unit_nan"), args_list)))
    } else if (descriptor %in% c("gaussian_diag", "continuous")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_diag"), args_list)))
    } else if (descriptor %in% c("gaussian_diag_nan", "continuous_nan")) {
      return(do.call(gaussian_model, c(list(type = "gaussian_diag_nan"), args_list)))

      # ------------------------------------------------------------------
      # Structural models - original names (retained for back-compatibility)
      # ------------------------------------------------------------------
    } else if (descriptor == "covariate") {
      return(do.call(covariate_model, args_list))
    } else if (descriptor == "distal_regression") {
      return(do.call(distal_regression_model, args_list))
    } else if (descriptor == "distal_pooled") {
      return(do.call(distal_pooled_model, args_list))
    } else if (descriptor == "distal_continuous") {
      return(do.call(distal_continuous_model, args_list))
    } else if (descriptor == "distal_continuous_regression") {
      return(do.call(distal_continuous_regression_model, args_list))
    } else if (descriptor == "distal_continuous_pooled") {
      return(do.call(distal_continuous_pooled_model, args_list))

      # ------------------------------------------------------------------
      # Structural models - user-friendly aliases
      #
      #   predict_class                -> covariate
      #   continuous_outcome           -> distal_continuous
      #   continuous_outcome_adjusted  -> distal_continuous_pooled
      #   continuous_outcome_moderated -> distal_continuous_regression
      #   categorical_outcome          -> distal_categorical
      #   categorical_outcome_adjusted -> distal_pooled
      #   categorical_outcome_moderated-> distal_regression
      # ------------------------------------------------------------------
    } else if (descriptor == "predict_class") {
      return(do.call(covariate_model, args_list))
    } else if (descriptor == "continuous_outcome") {
      return(do.call(distal_continuous_model, args_list))
    } else if (descriptor == "continuous_outcome_adjusted") {
      return(do.call(distal_continuous_pooled_model, args_list))
    } else if (descriptor == "continuous_outcome_moderated") {
      return(do.call(distal_continuous_regression_model, args_list))
    } else if (descriptor == "categorical_outcome") {
      return(do.call(distal_categorical_model, args_list))
    } else if (descriptor == "categorical_outcome_adjusted") {
      return(do.call(distal_pooled_model, args_list))
    } else if (descriptor == "categorical_outcome_moderated") {
      return(do.call(distal_regression_model, args_list))

    } else {
      stop(sprintf("Emission descriptor '%s' not recognized.", descriptor))
    }
  }
}
