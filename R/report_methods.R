# Methods-section reporting: report_methods() and grolts_report().
#
# Neither function computes anything new. Every number either one prints is
# read off the fitted object by .report_facts(), the single assembly helper
# both share, so the two reports cannot quote different values for the same
# fit. A number that would have to be computed here for the first time belongs
# in the fitting code instead, where print() and summary() can show it too.
#
# predict() and simulate() add nothing here: they are uses of a fit, and the
# fit does not record what is later done with it.

# The references the paragraphs may cite, keyed so that each paragraph names the
# ones it uses and the report prints only those, once, at the end.
.report_references <- c(
  dempster1977 = paste0(
    "Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood ",
    "from incomplete data via the EM algorithm. Journal of the Royal ",
    "Statistical Society: Series B, 39(1), 1-38."),
  galindo2006 = paste0(
    "Galindo Garre, F., & Vermunt, J. K. (2006). Avoiding boundary estimates ",
    "in latent class analysis by Bayesian posterior mode estimation. ",
    "Behaviormetrika, 33(1), 43-59. https://doi.org/10.2333/bhmk.33.43"),
  chen2008 = paste0(
    "Chen, J., Tan, X., & Zhang, R. (2008). Inference for normal mixtures in ",
    "mean and variance. Statistica Sinica, 18(2), 443-465."),
  tein2013 = paste0(
    "Tein, J.-Y., Coxe, S., & Cham, H. (2013). Statistical power to detect ",
    "the correct number of classes in latent profile analysis. Structural ",
    "Equation Modeling, 20(4), 640-657. ",
    "https://doi.org/10.1080/10705511.2013.824781"),
  masyn2013 = paste0(
    "Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. ",
    "In T. D. Little (Ed.), The Oxford Handbook of Quantitative Methods ",
    "(Vol. 2, pp. 551-611). Oxford University Press."),
  sorgente2025 = paste0(
    "Sorgente, A., Caliciuri, R., Robba, M., Lanz, M., & Zumbo, B. D. (2025). ",
    "A systematic review of latent class analysis in psychology: Examining ",
    "the gap between guidelines and research practice. Behavior Research ",
    "Methods, 57, 301."),
  vermunt2010 = paste0(
    "Vermunt, J. K. (2010). Latent class modeling with covariates: Two ",
    "improved three-step approaches. Political Analysis, 18(4), 450-469. ",
    "https://doi.org/10.1093/pan/mpq025"),
  bolck2004 = paste0(
    "Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent ",
    "structure models with categorical variables: One-step versus three-step ",
    "estimators. Political Analysis, 12(1), 3-27. ",
    "https://doi.org/10.1093/pan/mph001"),
  bakk2014 = paste0(
    "Bakk, Z., Oberski, D. L., & Vermunt, J. K. (2014). Relating latent ",
    "class assignments to external variables: Standard errors for correct ",
    "inference. Political Analysis, 22(4), 520-540. ",
    "https://doi.org/10.1093/pan/mpu003"),
  bakk2018 = paste0(
    "Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between ",
    "latent classes and external variables. Psychometrika, 83(4), 871-892. ",
    "https://doi.org/10.1007/s11336-017-9592-7"),
  vandeschoot2017 = paste0(
    "van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., & ",
    "Vermunt, J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on ",
    "latent trajectory studies. Structural Equation Modeling, 24(3), 451-467. ",
    "https://doi.org/10.1080/10705511.2016.1247646")
)

.fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
.pct  <- function(x) paste0(formatC(100 * x, format = "f", digits = 1), "%")
.plural <- function(n, one, many = paste0(one, "s")) if (n == 1) one else many

# The words a methods section uses for each indicator type, keyed on every
# spelling fit_mixture() accepts for it.
.indicator_word <- function(m) {
  m <- sub("_nan$", "", m)
  switch(m,
         binary = , bernoulli = "binary",
         categorical = , multinoulli = "polytomous",
         ordinal = "ordinal",
         count = , poisson = "count",
         continuous = , gaussian_diag = , gaussian = , gaussian_unit = ,
         structured_normal = "continuous",
         m)
}

# Which Bayes constant acts on which indicator type. A constant that has no
# parameter to act on is not mentioned: "every sentence must be true of this
# fit" rules out stating a prior that did nothing.
.prior_blocks <- function(types) {
  unique(c("latent", unlist(lapply(types, function(t) switch(t,
    binary = , binomial = , polytomous = , ordinal = "categorical",
    count = "poisson",
    continuous = "variances",
    NULL)))))
}

# The single place every reported fact is read off the object. `comparison` is
# the mixture_comparison the fit came from, when it came from one: a fit built
# inside compare_mixtures() carries no call of its own, and the call that
# rebuilds it is the comparison's.
.report_facts <- function(x, comparison = NULL) {
  lta    <- inherits(x, "lta_model")
  step1  <- if (lta && .is_threestep_lta(x)) x$step1 else x
  m      <- if (!is.null(x$step1_metrics)) x$step1_metrics else step1$metrics
  family <- if (inherits(x, "gmm")) "gmm" else if (inherits(x, "lcga")) "lcga"
            else if (inherits(x, "rmlca")) "rmlca" else if (lta) "lta"
            else "mixture"

  md       <- x$missing_data %||% list()
  n_input  <- md$n_input_rows %||% NROW(x$data)
  n_removed <- md$n_empty_rows %||% 0L
  dat      <- as.matrix(x$data)
  n_miss   <- md$n_missing %||% sum(is.na(dat))
  n_cells  <- md$n_cells %||% length(dat)
  w        <- x$sample_weights %||% x$weights_vec

  types <- if (family %in% c("gmm", "lcga")) {
    if (identical(x$growth$family, "gaussian")) "continuous"
    else .indicator_word(x$growth$family)
  } else if (!is.null(x$longitudinal$measurement)) {
    .indicator_word(x$longitudinal$measurement)
  } else {
    d <- x$measurement_descriptor
    if (is.character(d) && length(d) == 1L) .indicator_word(d)
    else tryCatch(unique(vapply(d, function(b) .indicator_word(
      if (is.list(b)) b$type %||% b[[1L]] else b), character(1))),
      error = function(e) "mixed-type")
  }

  shares <- if (lta) unname(x$prevalences[1L, ]) else unname(x$weights)
  n_eff  <- x$n_eff
  ave_pp <- separation <- shift <- NULL
  if (!lta && !is.null(x$log_resp) && length(shares) > 1L) {
    utils::capture.output(cd <- classification_diagnostics(x))
    ave_pp     <- diag(cd$ave_pp)
    separation <- attr(cd$separation, "min")
    shift      <- cd$class_shift
  }

  bc     <- x$bayes_constants %||% list()
  blocks <- .prior_blocks(types)
  bc     <- unlist(bc[intersect(blocks, names(bc))])

  # A fit made inside a comparison carries the comparison's internal call, in
  # its variable names; the line that rebuilds it is the comparison's own.
  call <- if (!is.null(comparison)) comparison$call else x$call
  seed <- if (!is.null(comparison)) comparison$call$random_state else x$random_state

  list(
    family      = family,
    unit        = if (lta) "status" else "class",
    k           = if (lta) x$n_statuses else x$n_components,
    n_input     = n_input,
    n_removed   = n_removed,
    n_analysed  = n_input - n_removed,
    n_eff       = n_eff,
    weighted    = !is.null(w) && any(w != 1) &&
                  !identical(x$weight_type, "frequency"),
    frequency   = identical(x$weight_type, "frequency"),
    survey      = isTRUE(x$has_survey_design),
    n_items     = if (!is.null(x$longitudinal$n_items)) x$longitudinal$n_items
                  else if (!is.null(x$growth))
                    length(x$growth$time_scores) else NCOL(x$data),
    n_times     = x$longitudinal$n_times %||% x$n_times,
    types       = types,
    n_steps     = x$n_steps %||% 1L,
    correction  = x$correction,
    assignment  = x$assignment,
    structural  = if (inherits(x$sm, "covariate")) "predictors"
                  else if (any(grepl("^distal", class(x$sm)))) "outcome"
                  else NULL,
    requested   = m$n_requested %||% m$n_starts,
    starts      = m$n_starts,
    replicated  = m$n_replicated,
    unreplicated = .is_unreplicated(m),
    refined_from = isTRUE(x$refined_from),
    converged   = x$converged,
    n_iter      = x$n_iter,
    priors      = bc,
    any_missing = n_miss > 0L,
    n_missing   = n_miss,
    n_cells     = n_cells,
    n_items_missing = md$n_items_affected %||% sum(colSums(is.na(dat)) > 0L),
    missing_rule = md$handled_by %||% "FIML (MAR assumption)",
    ll = m$ll, n_params = m$n_params, aic = m$aic, bic = m$bic,
    sabic = m$sabic, entropy = m$entropy,
    shares      = shares,
    n_class     = shares * n_eff,
    ave_pp      = ave_pp,
    separation  = separation,
    shift       = shift,
    growth      = x$growth,
    lta         = if (lta) list(
                    invariance = x$longitudinal$time_invariance,
                    homogeneous = isTRUE(x$tau_homogeneous),
                    mover_stayer = isTRUE(x$mover_stayer),
                    random_intercept = !is.null(x$ri),
                    n_classes = x$n_classes %||% 1L),
    version     = as.character(utils::packageVersion("mixtureEM")),
    r_version   = paste(R.version$major, R.version$minor, sep = "."),
    seed        = seed,
    call        = call,
    from_comparison = !is.null(comparison)
  )
}

# How many cases the ML third step moved, as one sentence appended to the
# structural paragraph; empty when the fit carries no step-3 posterior.
.shift_sentence <- function(sh) {
  if (is.null(sh)) return("")
  sprintf(paste0(" The modal class of %s %s (%s) changed between the ",
                 "measurement model and the third step."),
          .fmt(sh$n_moved, 0), .plural(round(sh$n_moved), "case", "cases"),
          .pct(sh$prop_moved))
}

# The model's name in the literature's words.
.report_model_name <- function(f) {
  g <- f$growth
  shape <- function() {
    deg <- g$degree %||% 1L
    c("linear", "quadratic", "cubic")[min(deg, 3L)]
  }
  # fit_gmm() with no random effect is the latent class growth model with its
  # residual structure free, and is named as one.
  fam <- if (f$family == "gmm" && identical(g$random_effects, "none")) "lcga"
         else f$family
  switch(fam,
    gmm   = sprintf("%d-class growth mixture model (%s growth, random %s)",
                    f$k, shape(), gsub("_", " and ", g$random_effects)),
    lcga  = sprintf("%d-class latent class growth analysis (%s growth, no within-class growth-factor variance)",
                    f$k, shape()),
    rmlca = sprintf("%d-class repeated-measures latent class model", f$k),
    lta   = sprintf("latent transition model with %d statuses over %d occasions",
                    f$k, f$n_times),
    {
      what <- if (all(f$types == "continuous")) "latent profile model"
              else if (any(f$types == "continuous")) "mixture model"
              else "latent class model"
      sprintf("%d-class %s", f$k, what)
    })
}

.para_sample <- function(f) {
  s <- sprintf("The analysis sample comprised %d cases", f$n_analysed)
  if (f$n_removed > 0L)
    s <- paste0(s, sprintf(" (of %d supplied; %d %s no observed indicator and %s removed)",
                           f$n_input, f$n_removed,
                           .plural(f$n_removed, "case had", "cases had"),
                           .plural(f$n_removed, "was", "were")))
  s <- paste0(s, ". ")
  ind <- if (f$family %in% c("gmm", "lcga"))
    sprintf("The outcome was measured on %d occasions (time scores %s).",
            f$n_items, paste(f$growth$time_scores, collapse = ", "))
  else if (f$family %in% c("lta", "rmlca"))
    sprintf("Each occasion was measured by %d %s indicators, over %d occasions.",
            f$n_items, paste(f$types, collapse = " and "), f$n_times)
  else
    sprintf("Classes were defined by %d %s indicators.", f$n_items,
            paste(f$types, collapse = " and "))
  s <- paste0(s, ind)
  if (f$frequency)
    s <- paste0(s, " Rows were weighted by frequency counts.")
  else if (f$weighted || f$survey)
    s <- paste0(s, " Estimation used sampling weights",
                if (f$survey) " and the survey design (strata and clusters)" else "",
                ".")
  list(text = s, refs = character())
}

.para_estimation <- function(f) {
  refs <- "dempster1977"
  s <- sprintf("A %s was estimated by maximum likelihood with the EM algorithm (Dempster, Laird, & Rubin, 1977)",
               .report_model_name(f))
  p <- f$priors
  if (length(p) && all(p == 0)) {
    s <- paste0(s, ", with no prior on any parameter. ")
  } else if (length(p)) {
    what <- c(latent = if (f$unit == "status")
                "initial-status and transition probabilities"
              else "class proportions",
              categorical = "item-response probabilities",
              poisson = "count rates",
              variances = "class-specific variances")
    on  <- p[p != 0]
    src <- c(if (any(names(on) != "variances")) "Galindo Garre & Vermunt, 2006",
             if ("variances" %in% names(on)) "Chen, Tan, & Zhang, 2008")
    refs <- c(refs, if (any(names(on) != "variances")) "galindo2006",
              if ("variances" %in% names(on)) "chen2008")
    parts <- sprintf("%s pseudo-%s for the %s", format(on),
                     ifelse(on == 1, "observation", "observations"),
                     what[names(on)])
    s <- paste0(s, sprintf(paste0(
      ", with weak priors that keep estimates off the boundary of the ",
      "parameter space (%s): %s"), paste(src, collapse = "; "),
      if (length(parts) > 1L)
        paste(paste(parts[-length(parts)], collapse = ", "), "and", parts[length(parts)])
      else parts))
    off <- names(p)[p == 0]
    s <- paste0(s, if (length(off)) sprintf("; none on the %s. ",
                                           paste(what[off], collapse = " or the "))
                   else ". ")
  } else s <- paste0(s, ". ")
  list(text = paste0(s, .starts_sentence(f)), refs = refs)
}

# The random starts, their replication and convergence: GRoLTS item 9 on its
# own, and the second half of the estimation paragraph.
.starts_sentence <- function(f) {
  s <- if (f$refined_from) {
    "The estimates were refined from a supplied starting solution. "
  } else if (!is.null(f$starts) && is.finite(f$starts) && f$starts >= 2L) {
    run <- if (!is.null(f$requested) && f$requested != f$starts)
      sprintf("The search requested %d random starts and carried %d to convergence",
              as.integer(f$requested), as.integer(f$starts))
    else sprintf("The search used %d random starts", as.integer(f$starts))
    rep <- if (!is.finite(f$replicated))
      "; none ran to convergence. "
    else if (f$replicated == 1L)
      "; the best log-likelihood was reached by one start only and is not replicated. "
    else sprintf("; the best log-likelihood was reached by %d of them. ",
                 as.integer(f$replicated))
    paste0(run, rep)
  } else "A single start was used. "
  paste0(s, if (isTRUE(f$converged))
    sprintf("The algorithm converged in %d iterations.", as.integer(f$n_iter))
  else sprintf("The algorithm had not converged when it stopped after %d iterations.",
               as.integer(f$n_iter)))
}

.para_missing <- function(f) {
  if (!f$any_missing) return(NULL)
  list(text = sprintf(paste0(
    "%d of %d indicator values (%s) were missing, in %d %s. Missing values ",
    "were handled by full-information maximum likelihood, which assumes they ",
    "are missing at random."),
    f$n_missing, f$n_cells, .pct(f$n_missing / f$n_cells),
    f$n_items_missing, .plural(f$n_items_missing,
      if (f$family %in% c("gmm", "lcga")) "occasion" else "indicator")),
    refs = character())
}

.para_fit <- function(f) {
  s <- sprintf(paste0("The model's log-likelihood was %s with %d free ",
                      "parameters (AIC = %s, BIC = %s, sample-size-adjusted BIC = %s)."),
               .fmt(f$ll), as.integer(f$n_params), .fmt(f$aic), .fmt(f$bic),
               .fmt(f$sabic))
  list(text = s, refs = character())
}

.para_classification <- function(f) {
  if (f$k < 2L) return(NULL)
  i <- which.min(f$shares)
  s <- sprintf("The smallest %s held %s of the sample%s (about %s cases).",
               f$unit, .pct(f$shares[i]),
               if (f$unit == "status") " at the first occasion" else "",
               .fmt(f$n_class[i], 0))
  if (!is.null(f$entropy) && is.finite(f$entropy))
    s <- paste0(sprintf("Relative entropy was %s. ", .fmt(f$entropy, 3)), s)
  if (!is.null(f$ave_pp))
    s <- paste0(s, sprintf(paste0(" The average posterior probability of the ",
                                  "assigned class ranged from %s to %s across classes."),
                           .fmt(min(f$ave_pp), 3), .fmt(max(f$ave_pp), 3)))
  s <- paste0(s, " These describe how cleanly cases are classified; they are ",
              "not criteria for the number of ", paste0(f$unit, "es"),
              " (Masyn, 2013).")
  refs <- "masyn2013"
  if (!is.null(f$separation)) {
    s <- paste0(s, sprintf(paste0(" The smallest Mahalanobis distance between ",
                                  "two class means, against the pooled ",
                                  "within-class variance, was %s; class ",
                                  "separation, more than sample size, governs ",
                                  "the power to detect the number of classes ",
                                  "(Tein, Coxe, & Cham, 2013)."),
                           .fmt(f$separation, 2)))
    refs <- c(refs, "tein2013")
  }
  list(text = s, refs = refs)
}

.para_structural <- function(f) {
  if (is.null(f$structural) && f$n_steps == 1L) return(NULL)
  what <- switch(f$structural %||% "", predictors = "class membership on the predictors",
                 outcome = "the distal outcome on class membership",
                 "the external variables on class membership")
  if (f$n_steps == 1L)
    return(list(text = sprintf(paste0(
      "The model relating %s was estimated jointly with the measurement ",
      "model (one-step)."), what), refs = character()))
  if (f$n_steps == 2L)
    return(list(text = sprintf(paste0(
      "The model relating %s was estimated by the two-step method: the ",
      "measurement model was held at its first-step estimates (Bakk & Kuha, 2018)."),
      what), refs = "bakk2018"))
  corr <- f$correction %||% "none"
  if (identical(corr, "BCH"))
    return(list(text = sprintf(paste0(
      "The model relating %s was estimated by the bias-adjusted three-step ",
      "method with the BCH correction (Bolck, Croon, & Hagenaars, 2004; ",
      "Vermunt, 2010), using %s assignment."), what, f$assignment),
      refs = c("bolck2004", "vermunt2010")))
  if (identical(corr, "ML"))
    return(list(text = sprintf(paste0(
      "The model relating %s was estimated by the bias-adjusted three-step ",
      "method with the maximum-likelihood correction (Vermunt, 2010), using ",
      "%s assignment, with standard errors that account for the first ",
      "step (Bakk, Oberski, & Vermunt, 2014).%s"), what, f$assignment,
      .shift_sentence(f$shift)),
      refs = c("vermunt2010", "bakk2014")))
  list(text = sprintf(paste0(
    "The model relating %s was estimated in three steps without a correction ",
    "for classification error, using %s assignment."), what, f$assignment),
    refs = character())
}

.para_lta <- function(f) {
  if (is.null(f$lta)) return(NULL)
  l <- f$lta
  s <- sprintf("Item-response probabilities were %s across occasions",
               switch(l$invariance %||% "full", full = "held equal",
                      none = "estimated freely", "held partly equal"))
  s <- paste0(s, if (l$homogeneous) ", and the transition probabilities were held equal over time"
                 else ", and the transition probabilities were estimated for each pair of occasions")
  extra <- c(if (l$mover_stayer) "a mover-stayer structure",
             if (l$random_intercept) "a random intercept",
             if (l$n_classes > 1L) sprintf("%d latent classes above the statuses", l$n_classes))
  if (length(extra)) s <- paste0(s, "; the model included ", paste(extra, collapse = " and "))
  list(text = paste0(s, "."), refs = character())
}

.para_repro <- function(f) {
  s <- sprintf("Models were estimated with the R package mixtureEM (version %s) in R %s.",
               f$version, f$r_version)
  s <- paste0(s, if (!is.null(f$seed))
    sprintf(" Random starts were drawn with random_state = %s.", format(f$seed))
  else " No random_state was set, so the random starts cannot be reproduced exactly.")
  cl <- if (!is.null(f$call)) paste(deparse(f$call, width.cutoff = 500L), collapse = " ")
  # A call built by do.call() carries its data inline and would print as
  # pages of numbers; it is the one case where no line is better than the line.
  if (!is.null(cl) && nchar(cl) <= 400L)
    s <- paste0(s, if (f$from_comparison) " The models are refitted by: "
                   else " The model is refitted by: ", cl)
  list(text = s, refs = character())
}

.para_enumeration <- function(cmp, k_reported) {
  tab  <- cmp$fit_table
  ks   <- tab$Classes
  unit <- if (identical(attr(cmp, "unit"), "statuses")) "statuses" else "classes"
  s <- sprintf("Models with %s %s were compared", paste(range(ks), collapse = " to "), unit)
  s <- paste0(s, if (1L %in% ks) ", including the one-class model. "
                 else "; the comparison did not include a one-class model. ")
  s <- paste0(s, "The comparison reported the log-likelihood, AIC, BIC, ",
              "sample-size-adjusted BIC, CAIC, AIC3, AWE and relative entropy for each. ")
  s <- paste0(s, sprintf("BIC was lowest for K = %d", cmp$best_k))
  s <- paste0(s, if (identical(as.integer(k_reported), as.integer(cmp$best_k)))
    ", the model reported here. "
  else sprintf("; the model reported here has K = %d. ", k_reported))
  unrep <- ks[which(tab$Unreplicated %in% TRUE)]
  if (length(unrep))
    s <- paste0(s, sprintf("The best log-likelihood was not replicated for K = %s. ",
                           paste(unrep, collapse = ", ")))
  s <- paste0(s, "Entropy is a classification diagnostic, not a criterion for ",
              "the number of ", unit, " (Masyn, 2013; Sorgente et al., 2025).")
  list(text = s, refs = c("masyn2013", "sorgente2025"))
}

.pick_model <- function(cmp, k) {
  k <- k %||% cmp$best_k
  fit <- cmp$models[[paste0("K", k)]]
  if (is.null(fit))
    stop(sprintf("The comparison holds no K = %d model; it holds K = %s.", k,
                 paste(cmp$fit_table$Classes, collapse = ", ")), call. = FALSE)
  list(fit = fit, k = k)
}

.methods_report <- function(paras) {
  paras <- Filter(Negate(is.null), paras)
  out   <- vapply(paras, `[[`, character(1), "text")
  keys  <- unique(unlist(lapply(paras, `[[`, "refs")))
  attr(out, "references") <- unname(sort(.report_references[keys]))
  class(out) <- "methods_report"
  out
}

#' Write the Methods Paragraph for a Fitted Mixture Model
#'
#' @description
#' Returns the paragraphs a methods section needs to describe a fitted model:
#' the sample, the estimator and its random starts, missing data, the
#' fit indices, classification quality, any structural model, and the line of
#' R that refits it. Every number is read off the fitted object, so the text is
#' true of the fit it was written from and of no other.
#'
#' @details
#' **Why this exists.** Sorgente et al. (2025) coded 313 latent class analyses
#' published in psychology between 2014 and 2024 against the published
#' tutorials: 23% reported the estimator, 16% the random starts, 17% how
#' missing data were handled, and 1.6% their analysis code. Every one of those
#' items is a value the package already stores; this function puts them in the
#' order and the wording a methods section uses.
#'
#' **What it does not do.** It reports and never grades. No sentence calls an
#' entropy, a class size or a replication count adequate: where the literature
#' attaches a cutoff, the cutoff is discussed in the help page that computes
#' the number (for entropy, [compare_mixtures()]). It does not say why the
#' number of classes is substantively right, what the classes mean, or whether
#' the sample suits the question; those are the author's.
#'
#' **Passing a comparison** (from [compare_mixtures()] or
#' [compare_longitudinal()]) adds the enumeration paragraph: the range of K,
#' whether the one-class model was in it, the indices reported, and which K had
#' the lowest BIC. The other paragraphs then describe the model with `k`
#' classes, by default the lowest-BIC one.
#'
#' Citations are returned as the `references` attribute and printed once at
#' the end, so the paragraphs serve both in-text and numbered citation styles.
#'
#' @param object A fitted model ([fit_mixture()], [fit_gmm()], [fit_lcga()],
#'   [fit_rmlca()], [fit_lta()]) or a comparison from [compare_mixtures()] or
#'   [compare_longitudinal()].
#' @param k For a comparison only: the number of classes of the model the
#'   paragraphs describe. Default is the comparison's lowest-BIC model.
#' @param ... Unused.
#'
#' @return A character vector of class `methods_report`, one paragraph per
#'   element, with the cited works in `attr(, "references")`. Printing it
#'   wraps the paragraphs and lists the references.
#'
#' @references
#' Sorgente, A., Caliciuri, R., Robba, M., Lanz, M., & Zumbo, B. D. (2025). A
#' systematic review of latent class analysis in psychology: Examining the gap
#' between guidelines and research practice. \emph{Behavior Research Methods},
#' \emph{57}, 301.
#'
#' Masyn, K. E. (2013). Latent class analysis and finite mixture modeling. In
#' T. D. Little (Ed.), \emph{The Oxford Handbook of Quantitative Methods}
#' (Vol. 2, pp. 551-611). Oxford University Press.
#'
#' Weller, B. E., Bowen, N. K., & Faubert, S. J. (2020). Latent class
#' analysis: A guide to best practice. \emph{Journal of Black Psychology},
#' \emph{46}(4), 287-311. \doi{10.1177/0095798420930932}
#'
#' Nylund-Gibson, K., & Choi, A. Y. (2018). Ten frequently asked questions
#' about latent class analysis. \emph{Translational Issues in Psychological
#' Science}, \emph{4}(4), 440-461. \doi{10.1037/tps0000176}
#'
#' @seealso [grolts_report()] for the trajectory-model checklist.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(600, 1, 0.5), nrow = 100)
#' fit <- fit_mixture(X, n_classes = 2, measurement = "binary", n_init = 5,
#'                    random_state = 1)
#' report_methods(fit)
#'
#' @export
report_methods <- function(object, ...) UseMethod("report_methods")

#' @rdname report_methods
#' @export
report_methods.mixture_model <- function(object, ...) {
  .report_methods_fit(.report_facts(object))
}

#' @rdname report_methods
#' @export
report_methods.lta_model <- function(object, ...) {
  .report_methods_fit(.report_facts(object))
}

#' @rdname report_methods
#' @export
report_methods.mixture_comparison <- function(object, k = NULL, ...) {
  pick <- .pick_model(object, k)
  f    <- .report_facts(pick$fit, comparison = object)
  .report_methods_fit(f, enumeration = .para_enumeration(object, pick$k))
}

.report_methods_fit <- function(f, enumeration = NULL) {
  .methods_report(list(.para_sample(f), .para_estimation(f), .para_missing(f),
                       enumeration, .para_fit(f), .para_lta(f),
                       .para_classification(f), .para_structural(f),
                       .para_repro(f)))
}

#' @export
print.methods_report <- function(x, ...) {
  for (p in x) cat(strwrap(p, width = 0.9 * getOption("width")), "", sep = "\n")
  refs <- attr(x, "references")
  if (length(refs)) {
    cat("References\n")
    for (r in refs) cat(strwrap(r, width = 0.9 * getOption("width"),
                                exdent = 2), sep = "\n")
  }
  invisible(x)
}

# The 21 GRoLTS items (van de Schoot et al., 2017, Table 1), worded as the
# checklist words them.
.grolts_items <- c(
  "1"   = "Is the metric of time used in the statistical model reported?",
  "2"   = "Is information presented about the mean and variance of time within a wave?",
  "3a"  = "Is the missing data mechanism reported?",
  "3b"  = "Is a description provided of what variables are related to attrition/missing data?",
  "3c"  = "Is a description provided of how missing data in the analyses were dealt with?",
  "4"   = "Is information about the distribution of the observed variables included?",
  "5"   = "Is the software mentioned?",
  "6a"  = "Are alternative specifications of within-class heterogeneity considered and clearly documented?",
  "6b"  = "Are alternative specifications of the between-class differences in variance-covariance matrix structure considered and clearly documented?",
  "7"   = "Are alternative shape/functional forms of the trajectories described?",
  "8"   = "If covariates have been used, can analyses still be replicated?",
  "9"   = "Is information reported about the number of random start values and final iterations included?",
  "10"  = "Are the model comparison (and selection) tools described from a statistical perspective?",
  "11"  = "Are the total number of fitted models reported, including a one-class solution?",
  "12"  = "Are the number of cases per class reported for each model?",
  "13"  = "If classification of cases in a trajectory is the goal, is entropy reported?",
  "14a" = "Is a plot included with the estimated mean trajectories of the final solution?",
  "14b" = "Are plots included with the estimated mean trajectories for each model?",
  "14c" = "Is a plot included of the estimated means of the final model with the observed individual trajectories, by class?",
  "15"  = "Are characteristics of the final class solution numerically described?",
  "16"  = "Are the syntax files available?"
)

#' Answer the GRoLTS Reporting Checklist for a Trajectory Model
#'
#' @description
#' Walks the 21 items of the GRoLTS checklist (Guidelines for Reporting on
#' Latent Trajectory Studies; van de Schoot et al., 2017) for a fitted growth
#' model or a comparison of them, and sorts each item into one of three
#' groups: answered from the fit, which you must supply, and not applicable to
#' this model.
#'
#' @details
#' Van de Schoot et al. (2017) applied the checklist to 38 published
#' trajectory studies; only one reported the number of random starts (item
#' 9), and none showed the mean trajectories of every model compared (item
#' 14b). Most of what the checklist asks for is already on the fitted object;
#' this function gathers it, and names what the fit cannot answer: the
#' missing-data mechanism, the variables related to attrition, the
#' distribution of the observed variables, the plots, and which alternative
#' specifications were considered and why.
#'
#' It deliberately prints no score. A count of items met invites treating the
#' items as interchangeable, which they are not.
#'
#' A comparison from [compare_longitudinal()] answers the enumeration items
#' (10, 11, 12) for every model; a single fit answers them for itself only, and
#' says so.
#'
#' @param object A [fit_gmm()] or [fit_lcga()] fit, or a comparison of them
#'   from [compare_longitudinal()].
#' @param k For a comparison only: the number of classes of the final model.
#'   Default is the comparison's lowest-BIC model.
#' @param ... Unused.
#'
#' @return A data frame of class `grolts_report` with one row per item and
#'   columns `item`, `question`, `status` (`"from the fit"`, `"you supply"` or
#'   `"not applicable"`) and `answer`. Printing it groups the items by status.
#'
#' @references
#' van de Schoot, R., Sijbrandij, M., Winter, S. D., Depaoli, S., & Vermunt,
#' J. K. (2017). The GRoLTS-checklist: Guidelines for reporting on latent
#' trajectory studies. \emph{Structural Equation Modeling}, \emph{24}(3),
#' 451-467. \doi{10.1080/10705511.2016.1247646}
#'
#' @seealso [report_methods()] for the methods paragraph.
#'
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(400), 100) + rep(c(0, 1, 2, 3), each = 100)
#' fit <- fit_lcga(y, times = 4, n_classes = 2, family = "gaussian",
#'                 n_init = 3, random_state = 1)
#' grolts_report(fit)
#'
#' @export
grolts_report <- function(object, ...) UseMethod("grolts_report")

#' @rdname grolts_report
#' @export
grolts_report.mixture_model <- function(object, ...) {
  if (!inherits(object, c("gmm", "lcga")))
    stop("grolts_report() is for trajectory models: a fit_gmm() or ",
         "fit_lcga() fit, or compare_longitudinal() over one. For any other ",
         "model use report_methods().", call. = FALSE)
  .grolts(object, NULL)
}

#' @rdname grolts_report
#' @export
grolts_report.mixture_comparison <- function(object, k = NULL, ...) {
  pick <- .pick_model(object, k)
  if (!inherits(pick$fit, c("gmm", "lcga")))
    stop("grolts_report() is for trajectory models; this comparison holds ",
         class(pick$fit)[1L], " fits. Use report_methods().", call. = FALSE)
  .grolts(pick$fit, object)
}

.grolts <- function(fit, cmp) {
  f <- .report_facts(fit, comparison = cmp)
  g <- f$growth
  gmm <- f$family == "gmm" && !identical(g$random_effects, "none")
  FIT <- "from the fit"; YOU <- "you supply"; NA_ <- "not applicable"
  sizes <- function(x) paste(.pct(x), collapse = ", ")
  starts <- .starts_sentence(f)

  ans <- list(
    "1"  = list(FIT, sprintf("Time scores %s for occasions %s.",
                             paste(g$time_scores, collapse = ", "),
                             paste(g$time_labels, collapse = ", "))),
    "2"  = list(YOU, paste0("The model uses one fixed time score per occasion; ",
                            "report how much actual measurement times varied within each wave.")),
    "3a" = list(YOU, paste0("The estimator assumes values are missing at random; ",
                            "argue whether that is plausible. mcar_test() tests the stronger MCAR assumption.")),
    "3b" = list(YOU, "Name the variables related to attrition or missingness."),
    "3c" = list(FIT, paste0((.para_missing(f) %||% list(text = sprintf(
                  "No values were missing on the %d occasions.", f$n_items)))$text,
                  if (f$n_removed > 0L) sprintf(
                    " %d of %d cases had no observed occasion and were removed.",
                    f$n_removed, f$n_input) else "")),
    "4"  = list(YOU, sprintf("Describe the observed %s outcome at each occasion (means, SDs, skew).",
                             f$types)),
    "5"  = list(FIT, sprintf("R package mixtureEM %s, R %s.", f$version, f$r_version)),
    "6a" = list(YOU, sprintf(paste0("This fit is a %s. Report which other within-class ",
                                    "specifications (e.g., LCGA against GMM) were fitted and why this one was kept."),
                             .report_model_name(f))),
    "6b" = if (gmm) list(YOU, sprintf(paste0(
                  "This fit holds the growth-factor covariance %s and the residual variances %s. ",
                  "Report which other structures were considered."),
                  if (isTRUE(g$psi_equal)) "equal across classes" else "free across classes",
                  if (isTRUE(g$residual_equal)) "equal across classes" else "free across classes"))
           else if (f$types == "continuous") list(YOU, sprintf(paste0(
                  "This fit has no growth-factor covariance; its residual variances are %s. ",
                  "Report which other structures were considered."),
                  if (isTRUE(g$residual_equal)) "equal across classes" else "free across classes"))
           else list(NA_, paste0("This model has no growth-factor covariance and no ",
                                 "residual variance to structure.")),
    "7"  = list(YOU, sprintf("This fit uses %s growth. Report the other shapes fitted, if any.",
                             c("linear", "quadratic", "cubic")[min(g$degree %||% 1L, 3L)])),
    "8"  = if (length(g$covariate_names)) list(FIT, sprintf(
                  "Covariates %s, with effects %s across classes.",
                  paste(g$covariate_names, collapse = ", "),
                  if (isTRUE(g$covariate_equal)) "equal" else "free"))
           else list(NA_, "No covariates were used."),
    "9"  = list(FIT, starts),
    "10" = if (!is.null(cmp)) list(FIT, .para_enumeration(cmp, f$k)$text)
           else list(YOU, "Pass the compare_longitudinal() result to answer this from the fits."),
    "11" = if (!is.null(cmp)) list(FIT, sprintf(
                  "%d models were fitted, K = %s; the one-class model %s.",
                  nrow(cmp$fit_table), paste(cmp$fit_table$Classes, collapse = ", "),
                  if (1L %in% cmp$fit_table$Classes) "was included" else "was not included"))
           else list(YOU, "Pass the compare_longitudinal() result to answer this from the fits."),
    "12" = if (!is.null(cmp)) list(FIT, paste(vapply(cmp$models, function(m)
                  sprintf("K = %d: %s", m$n_components, sizes(m$weights)), character(1)),
                  collapse = "; "))
           else list(FIT, sprintf("Final model only: %s. Pass the comparison for every model.",
                                  sizes(f$shares))),
    "13" = if (f$k > 1L) list(FIT, sprintf("Relative entropy %s.", .fmt(f$entropy, 3)))
           else list(NA_, "A one-class model classifies no one."),
    "14a" = list(YOU, "plot(fit) draws the estimated mean trajectories; include it."),
    "14b" = list(YOU, "Plot each compared model's mean trajectories, not only the final one's."),
    "14c" = list(YOU, "Overlay the observed individual trajectories on each class's estimated mean."),
    "15" = list(FIT, sprintf("Growth-factor means by class are in summary(fit); class shares %s.",
                             sizes(f$shares))),
    "16" = list(FIT, .para_repro(f)$text)
  )
  out <- data.frame(item = names(.grolts_items), question = unname(.grolts_items),
                    status = vapply(ans, `[[`, character(1), 1L),
                    answer = vapply(ans, `[[`, character(1), 2L),
                    stringsAsFactors = FALSE, row.names = NULL)
  class(out) <- c("grolts_report", "data.frame")
  out
}

#' @export
print.grolts_report <- function(x, ...) {
  w <- 0.9 * getOption("width")
  cat("GRoLTS checklist (van de Schoot et al., 2017)\n")
  for (st in c("from the fit", "you supply", "not applicable")) {
    rows <- x[x$status == st, , drop = FALSE]
    if (!nrow(rows)) next
    cat("\n", switch(st, "from the fit" = "Answered from the fit",
                     "you supply" = "You must supply",
                     "not applicable" = "Not applicable to this model"), "\n", sep = "")
    for (i in seq_len(nrow(rows))) {
      cat(strwrap(sprintf("%s. %s", rows$item[i], rows$question[i]), width = w,
                  exdent = 4), sep = "\n")
      cat(strwrap(rows$answer[i], width = w, indent = 4, exdent = 4), "", sep = "\n")
    }
  }
  invisible(x)
}
