# ==============================================================================
# S3 Emission Factory (List Dispatcher)
# ==============================================================================

# Call a model constructor, dropping any argument it cannot accept.
#
# Extra arguments reach build_emission() through fit_mixture()'s `...`, which is
# shared by the measurement and the structural model. A longitudinal measurement
# model needs n_items / n_times / sub_model, none of which a covariate or distal
# constructor accepts, so unusable arguments are dropped here rather than
# surfacing as an "unused argument" error from the wrong constructor.
.construct_emission <- function(constructor, args) {
  formal_names <- names(formals(constructor))
  if (!"..." %in% formal_names)
    args <- args[names(args) %in% formal_names]
  do.call(constructor, args)
}

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
      return(.construct_emission(categorical_model, c(list(type = "bernoulli"), args_list)))
    } else if (descriptor %in% c("bernoulli_nan", "binary_nan")) {
      return(.construct_emission(categorical_model, c(list(type = "bernoulli_nan"), args_list)))
    } else if (descriptor %in% c("multinoulli", "categorical")) {
      return(.construct_emission(categorical_model, c(list(type = "multinoulli"), args_list)))
    } else if (descriptor %in% c("multinoulli_nan", "categorical_nan")) {
      return(.construct_emission(categorical_model, c(list(type = "multinoulli_nan"), args_list)))
    } else if (descriptor %in% c("gaussian_unit", "gaussian")) {
      return(.construct_emission(gaussian_model, c(list(type = "gaussian_unit"), args_list)))
    } else if (descriptor %in% c("gaussian_unit_nan", "gaussian_nan")) {
      return(.construct_emission(gaussian_model, c(list(type = "gaussian_unit_nan"), args_list)))
    } else if (descriptor %in% c("gaussian_diag", "continuous")) {
      return(.construct_emission(gaussian_model, c(list(type = "gaussian_diag"), args_list)))
    } else if (descriptor %in% c("gaussian_diag_nan", "continuous_nan")) {
      return(.construct_emission(gaussian_model, c(list(type = "gaussian_diag_nan"), args_list)))
    } else if (descriptor %in% c("poisson", "count")) {
      return(.construct_emission(poisson_model, c(list(type = "poisson"), args_list)))
    } else if (descriptor %in% c("poisson_nan", "count_nan")) {
      return(.construct_emission(poisson_model, c(list(type = "poisson_nan"), args_list)))

      # ------------------------------------------------------------------
      # Longitudinal measurement model: J items observed at T occasions,
      # optionally with response parameters held equal across occasions.
      # Extra arguments (n_items, n_times, sub_model, invariant_items,
      # max_val) reach here through fit_mixture_internal()'s `...`.
      # ------------------------------------------------------------------
    } else if (descriptor == "time_blocks") {
      return(.construct_emission(time_blocks_model, args_list))

      # ------------------------------------------------------------------
      # Multiple-group measurement model: J items observed once per case,
      # laid out as G group-blocks with a case's own group's block populated
      # and every other block structurally missing. Extra arguments
      # (n_items, n_groups, sub_model, invariant_items, max_val) reach here
      # through fit_mixture_internal()'s `...`. See R/group_blocks.R.
      # ------------------------------------------------------------------
    } else if (descriptor == "group_blocks") {
      return(.construct_emission(group_blocks_model, args_list))

      # ------------------------------------------------------------------
      # Latent class growth model: one outcome at T occasions, each class
      # following its own polynomial in time on the link scale. Extra
      # arguments (design, family) reach here through fit_lcga(). See
      # R/lcga.R.
      # ------------------------------------------------------------------
    } else if (descriptor == "lcga") {
      return(.construct_emission(lcga_model, args_list))

      # ------------------------------------------------------------------
      # Growth mixture model: one continuous outcome at T occasions, each
      # class following its own polynomial in time, with cases varying about
      # their class's curve through random growth factors. Extra arguments
      # (design, random_effects, psi, residual, residual_equal) reach here
      # through fit_gmm(). See R/structured_normal.R.
      # ------------------------------------------------------------------
    } else if (descriptor %in% c("structured_normal", "gmm")) {
      return(.construct_emission(structured_normal_model, args_list))

      # ------------------------------------------------------------------
      # Structural models - original names (retained for back-compatibility)
      # ------------------------------------------------------------------
    } else if (descriptor == "covariate") {
      return(.construct_emission(covariate_model, args_list))
    } else if (descriptor == "distal_regression") {
      return(.construct_emission(distal_regression_model, args_list))
    } else if (descriptor == "distal_pooled") {
      return(.construct_emission(distal_pooled_model, args_list))
    } else if (descriptor == "distal_continuous") {
      return(.construct_emission(distal_continuous_model, args_list))
    } else if (descriptor == "distal_continuous_regression") {
      return(.construct_emission(distal_continuous_regression_model, args_list))
    } else if (descriptor == "distal_continuous_pooled") {
      return(.construct_emission(distal_continuous_pooled_model, args_list))

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
      return(.construct_emission(covariate_model, args_list))
    } else if (descriptor == "continuous_outcome") {
      return(.construct_emission(distal_continuous_model, args_list))
    } else if (descriptor == "continuous_outcome_adjusted") {
      return(.construct_emission(distal_continuous_pooled_model, args_list))
    } else if (descriptor == "continuous_outcome_moderated") {
      return(.construct_emission(distal_continuous_regression_model, args_list))
    } else if (descriptor == "categorical_outcome") {
      return(.construct_emission(distal_categorical_model, args_list))
    } else if (descriptor == "categorical_outcome_adjusted") {
      return(.construct_emission(distal_pooled_model, args_list))
    } else if (descriptor == "categorical_outcome_moderated") {
      return(.construct_emission(distal_regression_model, args_list))

    } else {
      stop(sprintf("Emission descriptor '%s' not recognized.", descriptor))
    }
  }
}
