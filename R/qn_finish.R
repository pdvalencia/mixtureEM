# ==============================================================================
# Newton-type finish after EM
# ==============================================================================
#
# EM climbs quickly while the classes are sorting themselves out and then slows
# to a crawl, and it crawls hardest exactly where categorical models most often
# end up: an optimum with some response probabilities at 0 or 1. The logits of
# those cells head to infinity, the likelihood is nearly flat along them, and
# the fraction of missing information that governs EM's rate approaches one in
# the directions coupled to them. On a four-class binary LCA with a class
# membership regression (N = 1373) the EM increment was still about 2e-6 per
# iteration two thousand iterations past the point where the classes had
# settled, and EM needed 45,000 iterations to meet an absolute rule of 1e-12.
#
# A stopping rule on the change in the log-likelihood cannot tell that crawl
# from convergence. The package's rule stops once a step gains less than
# `abs = 1e-6` or a relative 1e-8, which at |LL| ~ 4000 is ~4e-5 a step, and on
# the benchmark fits that left the reported optimum 0.002-0.019 below the
# maximum: enough to move a classification table by 2e-3 and a three-step
# estimate built on it by a few hundredths (internal/RECORDS.md, "Ridge
# convergence diagnosis").
#
# The remedy is the standard one: let EM find the basin, then finish with a
# Newton-type method, whose rate does not depend on the missing information,
# and judge convergence on the gradient rather than on the size of the last
# step. Here that is BHHH (Berndt, Hall, Hall and Hausman, 1974): Newton steps
# with the curvature taken as the outer product of the case-level scores, on
# the model's own packed parameter vector (the unconstrained packing
# R/step3_variance.R and R/lta_core.R already use for standard errors), on the
# same penalised objective EM climbs.
#
# The two details that decide whether the finish lands where EM would have:
#
#   * The curvature must be the data's. Near these optima some directions carry
#     almost no information at all, and a small interior maximum can sit on
#     one of them: a class-by-item logit at -9 with a standard error of 249,
#     worth 1e-5 of log-likelihood over the boundary. BFGS starts from an
#     identity curvature, and its line search carried such a coordinate from
#     -7 to -107 in one step, past the maximum and onto the boundary, where the
#     gradient vanishes and nothing brings it back; long EM runs, and Newton on
#     the observed information, both stop at -9. The outer product of scores is
#     the empirical information, costs one gradient's worth of likelihood
#     evaluations, and makes the step along such a direction as small as the
#     data say it should be.
#   * Steps are also capped at one unit in any coordinate, and a coordinate
#     that crosses |15| is taken to be heading for the boundary: it is sent to
#     |25| (a probability of 1e-11, the box R/lta_core.R's polish uses) and held
#     there. At |25| its remaining contribution is below 1e-8 of log-likelihood
#     at any sample size fitted here, and a free coordinate there only makes
#     the curvature singular.
#
# Measured from each fit's own EM winner, the finish reaches the fully
# converged maximum -- the one long EM runs reach -- on every benchmark fit, in
# seconds rather than the tens of thousands of iterations EM would need
# (RECORDS.md, same entry).

# Past this, on the logit scale, a coordinate is sent to `.qn_edge` and held.
.qn_wall <- 15
.qn_edge <- 25

# The largest change any coordinate may take in one step.
.qn_radius <- 1

# The iterations stop when every free gradient component is below `gtol_stop`
# per unit of case weight, when every step along the search direction would
# lower the objective, or at `maxit`. The rule is on the gradient, not on the gain:
# BHHH converges linearly, and near the maximum its steps gain less than the
# objective's own rounding (1e-12 at |LL| = 1e4, often exactly 0) while the
# parameters are still moving. A rule on the gain at 1e-15 of the
# objective stopped a three-step LTA with its estimates 1e-7 from the maximum,
# and one iteration's difference in where it stopped moved them by that much
# between a data set and the same data with every case duplicated; on the
# gradient the two agree to 1e-11. Relative to the case weight, so that
# duplication, which doubles the gradient, stops at the same point (the step
# is already invariant to it, the outer product doubling with the gradient).
# With an analytic gradient the rule is 1e-12; a central-difference gradient
# is only good to ~1e-7 at |LL| ~ 4000, so it gets 1e-9 (1.4e-6 at N = 1373).
.qn_gtol_stop    <- 1e-12
.qn_gtol_stop_fd <- 1e-9

# A fit whose finish ends with every free gradient component below this, per
# unit of case weight, is reported as converged, whether or not EM met its own
# rule first: EM's iteration cap is a statement about EM, not about where the
# fit ended up.
.qn_gtol <- 1e-6

.qn_fd_grad <- function(f, par) {
  h <- 1e-5 * pmax(1, abs(par))
  vapply(seq_along(par), function(i) {
    a <- par; a[i] <- a[i] + h[i]
    b <- par; b[i] <- b[i] - h[i]
    (f(a) - f(b)) / (2 * h[i])
  }, numeric(1))
}

# Case-level scores of `case_ll` by central differences, for the coordinates
# `idx` only: an n x length(idx) matrix.
.qn_fd_scores <- function(case_ll, par, idx) {
  h <- 1e-5 * pmax(1, abs(par[idx]))
  S <- vapply(seq_along(idx), function(m) {
    a <- par; a[idx[m]] <- a[idx[m]] + h[m]
    b <- par; b[idx[m]] <- b[idx[m]] - h[m]
    (case_ll(a) - case_ll(b)) / (2 * h[m])
  }, numeric(length(case_ll(par))))
  matrix(S, ncol = length(idx))
}

# The finish. `case_ll(par)` is the per-row log-likelihood, `w` the row
# weights, `prior(par)` the log-prior (0 when there is none), and
# `scores(par, idx)` the unweighted per-row scores for the coordinates `idx`
# (central differences of `case_ll` when not supplied). Returns NULL when
# nothing is free or the start is not finite, so every caller can fall back to
# the EM solution.
.qn_finish <- function(case_ll, w, par0, prior = function(p) 0,
                       scores = NULL, gtol_stop = .qn_gtol_stop,
                       maxit = 500L) {
  scores <- scores %||% function(p, idx) .qn_fd_scores(case_ll, p, idx)
  value  <- function(p) {
    v <- sum(w * case_ll(p)) + prior(p)
    if (is.finite(v)) v else -Inf
  }
  x  <- par0
  fx <- value(x)
  if (!is.finite(fx)) return(NULL)
  sw <- sum(w)
  free <- which(is.finite(x) & abs(x) <= .qn_wall)
  if (!length(free)) return(NULL)

  # Coordinates already past the wall go to the edge, if that does not cost.
  to_edge <- function(x, fx, idx) {
    e <- x
    e[idx] <- sign(x[idx]) * .qn_edge
    fe <- value(e)
    if (fe >= fx) list(x = e, fx = fe) else list(x = x, fx = fx)
  }
  held <- setdiff(seq_along(x), free)
  if (length(held)) { r <- to_edge(x, fx, held); x <- r$x; fx <- r$fx }

  it <- 0L
  g  <- NULL
  repeat {
    it <- it + 1L
    S  <- scores(x, free)
    fr <- function(v) { p <- x; p[free] <- v; prior(p) }
    gp <- .qn_fd_grad(fr, x[free])
    g  <- colSums(S * w) + gp
    if (max(abs(g)) < gtol_stop * sw || it > maxit) break

    B  <- crossprod(S * sqrt(w))
    e  <- eigen(B, symmetric = TRUE)
    lam <- pmax(e$values, 1e-10 * max(e$values, 1e-300))
    d  <- as.vector(e$vectors %*% (crossprod(e$vectors, g) / lam))
    big <- max(abs(d))
    if (big > .qn_radius) d <- d * .qn_radius / big

    a <- 1
    repeat {
      xn <- x
      xn[free] <- x[free] + a * d
      fn <- value(xn)
      if (fn >= fx) break
      a <- a / 2
      if (a < 1e-10) break
    }
    if (a < 1e-10) break

    hit <- free[abs(xn[free]) > .qn_wall]
    if (length(hit)) {
      r <- to_edge(xn, fn, hit); xn <- r$x; fn <- r$fx
      free <- setdiff(free, hit)
      if (!length(free)) { x <- xn; fx <- fn; break }
    }
    x  <- xn
    fx <- fn
  }
  list(par = x, value = fx,
       max_gradient = if (length(free)) max(abs(g)) else 0,
       converged = !length(free) || max(abs(g)) < .qn_gtol * sw,
       iterations = it, n_held = length(x) - length(free))
}

# Which of several EM solutions to finish. The restart that EM ranked first is
# not always the one with the highest maximum: restarts are stopped at
# different points of their crawl, so two basins within a few hundredths of
# each other can swap order once both are finished. Candidates within `band`
# of the best are finished in rank order, at most `max_finish` of them, and one
# is skipped when it is the same solution as one already taken. "Same" is read
# off the per-case log-likelihoods, which do not depend on how the classes are
# labelled: on the benchmark that calibrated this, restarts converging to one
# maximum differed by at most 0.006 per case before finishing, and the nearest
# distinct maximum by 0.02 (RECORDS.md, "Ridge convergence diagnosis").
.qn_pick <- function(scores, signatures, weights, band = 0.05,
                     max_finish = 3L, same = 0.01) {
  ord  <- order(scores, decreasing = TRUE)
  ord  <- ord[scores[ord] >= scores[ord[1L]] - band]
  keep <- integer(0)
  sw   <- sum(weights)
  for (i in ord) {
    if (length(keep) >= max_finish) break
    dup <- any(vapply(keep, function(j)
      sum(weights * abs(signatures[[i]] - signatures[[j]])) / sw < same,
      logical(1)))
    if (!dup) keep <- c(keep, i)
  }
  keep
}

# ------------------------------------------------------------------------------
# fit_mixture()
# ------------------------------------------------------------------------------

# The measurement families the finish is run for: categorical ones, where the
# boundary crawl lives and the objective is bounded. A Gaussian or count model
# is left to EM, whose rule was measured adequate there (RECORDS.md, "R3's
# W1"), and where an unpenalised variance has no bound to stop a line search.
.qn_mixture_families <- c("bernoulli", "multinoulli")

.qn_mixture_supported <- function(model_state, Y) {
  if (length(model_state$frozen) || model_state$n_components < 2L) return(FALSE)
  mm <- model_state$mm
  subs <- if (inherits(mm, c("blocks", "nested"))) mm$models else list(mm)
  ok <- vapply(subs, function(m)
    inherits(m, "bernoulli_dif") || .step1_family(m) %in% .qn_mixture_families,
    logical(1))
  if (!all(ok)) return(FALSE)
  # With structural data the structural model has to be the one the packing
  # covers: a class-membership regression. Anything else is left to EM.
  is.null(Y) || inherits(model_state$sm, "covariate")
}

# The Dirichlet prior the covariate M-step writes as pseudo-rows (.fit_mnl(),
# R/covariate.R): alpha / (K * U0) of a case in every class at every distinct
# covariate pattern. .em_log_prior() leaves the class probabilities to the
# structural model when it supplies them, so the term is priced here.
.qn_covariate_prior <- function(sm, Y) {
  alpha <- .bayes_alpha(sm, "latent")
  if (!(alpha > 0)) return(function(sm) 0)
  Z  <- complete_covariates(as.matrix(Y))
  Zm <- unique(if (isTRUE(sm$intercept)) cbind(1, Z) else Z)
  K  <- sm$n_components
  wt <- alpha / (K * nrow(Zm))
  function(sm) wt * sum(log(pmax(softmax_rows(Zm %*% t(sm$parameters$beta)),
                                 1e-300)))
}

.qn_finish_mixture <- function(model_state, X, Y) {
  use_sm <- !is.null(Y) && .supplies_class_probs(model_state$sm)
  pm <- .step1_pack_mm(model_state$mm)
  if (is.null(pm)) return(NULL)
  par0 <- if (use_sm) {
    ps <- .step1_pack_sm(model_state$sm)
    if (is.null(ps)) return(NULL)
    c(pm, ps)
  } else .step1_pack(model_state)
  nm <- length(pm)

  unpack <- function(par) {
    if (!use_sm) return(.step1_unpack(model_state, par))
    st <- model_state
    st$mm <- .step1_unpack_mm(st$mm, par[seq_len(nm)])
    st$sm <- .step1_unpack_sm(st$sm, par[-seq_len(nm)])
    st
  }
  Yp     <- if (use_sm) Y else NULL
  marg   <- tryCatch(.em_prior_marginals(model_state, X), error = function(e) NULL)
  cprior <- if (use_sm) .qn_covariate_prior(model_state$sm, Y) else function(sm) 0
  case_ll <- function(par) {
    st <- unpack(par)
    if (use_sm)
      logsumexp(log_likelihood(st$mm, X) + log_likelihood(st$sm, Y), MARGIN = 1)
    else
      logsumexp(sweep(log_likelihood(st$mm, X), 2,
                      log(pmax(st$weights, 1e-300)), "+"), MARGIN = 1)
  }
  prior <- function(par) {
    st <- unpack(par)
    lp <- .em_log_prior(st, X, Yp, marg)
    if (is.na(lp)) NA_real_ else lp + cprior(st$sm)
  }
  if (is.na(prior(par0))) return(NULL)
  w <- model_state$sample_weights
  start <- sum(w * case_ll(par0)) + prior(par0)
  if (!is.finite(start)) return(NULL)
  res <- .qn_finish(case_ll, w, par0, prior, gtol_stop = .qn_gtol_stop_fd)
  if (is.null(res) || res$value <= start) return(NULL)

  st <- unpack(res$par)
  e  <- e_step(st, X, Yp)
  st$log_resp    <- e$log_resp
  st$lower_bound <- e$log_prob_norm
  st$converged   <- isTRUE(model_state$converged) || res$converged
  st$qn_finish <- list(iterations = res$iterations, n_held = res$n_held,
                       max_gradient = res$max_gradient, converged = res$converged,
                       gain = res$value - start)
  st
}

# Finish the leading EM solutions of a multi-start search and return the best.
# `states` are the fitted restarts, `score` the objective the search ranked
# them on (fit_em()'s ll_of()).
.qn_finish_search <- function(states, score, X, Y) {
  states <- Filter(Negate(is.null), states)
  if (!length(states) || !.qn_mixture_supported(states[[1]], Y)) return(NULL)
  scores <- vapply(states, score, numeric(1))
  pick <- .qn_pick(scores, lapply(states, `[[`, "lower_bound"),
                   states[[1]]$sample_weights)
  best <- NULL; best_score <- -Inf
  for (i in pick) {
    st <- .qn_finish_mixture(states[[i]], X, Y) %||% states[[i]]
    s  <- score(st)
    if (s > best_score) { best <- st; best_score <- s }
  }
  best
}

# ------------------------------------------------------------------------------
# fit_lta()
# ------------------------------------------------------------------------------

# A latent transition model with categorical indicators and no free random
# intercept loading. A continuous random intercept has its own search and
# polish (R/lta.R), validated on their own benchmarks, and is left to them.
.qn_lta_supported <- function(state) {
  if (length(state$frozen) || isTRUE(state$mm_fixed)) return(FALSE)
  if ((state$n_classes %||% 1L) > 1L) return(FALSE)
  if (!is.null(state$ri) && .lta_ri_loading_free(state)) return(FALSE)
  if (!.lta_scores_supported(state, shared_ok = TRUE)) return(FALSE)
  class(state$mm$models[[1]])[1] %in%
    c("bernoulli", "bernoulli_nan", "ordinal", "ordinal_nan")
}

# The packed vector of an LTA state, reduced to the parameters it actually
# has. The transition-free model (`tau_independent`) ties every origin row of
# an occasion's matrix to one vector; with a plain `tau` the ordinary packing
# describes those rows one by one, so they are packed with the flag off and
# mapped back onto one free vector per occasion -- the reduction
# .lta_threestep_vcov() uses for the same model. With a transition regression
# the origin-free design already has one row per occasion and packs as it is.
.qn_lta_reduction <- function(state) {
  shared <- isTRUE(state$tau_independent) && is.null(state$tau_beta)
  st <- state
  if (shared) st$tau_independent <- FALSE
  lay  <- .lta_par_layout(st)
  full <- .lta_par_pack(st, lay)
  len  <- vapply(lay, function(b) as.integer(b$len), integer(1))
  if (sum(len) != length(full)) return(NULL)
  beg  <- cumsum(len) - len
  src  <- integer(length(full))
  n_red <- 0L
  first <- list()
  for (i in seq_along(lay)) {
    b    <- lay[[i]]
    cols <- beg[i] + seq_len(len[i])
    key  <- if (shared && identical(b$kind, "tau")) paste0("tau", b$i_mat)
    if (!is.null(key) && !is.null(first[[key]])) {
      src[cols] <- first[[key]]
    } else {
      src[cols] <- n_red + seq_len(len[i])
      n_red <- n_red + len[i]
      if (!is.null(key)) first[[key]] <- src[cols]
    }
  }
  r0 <- full[match(seq_len(n_red), src)]
  if (max(abs(r0[src] - full)) > 1e-8) return(NULL)
  list(state = st, layout = lay, src = src, r0 = r0, shared = shared)
}

.qn_finish_lta <- function(state, X, alpha, marginals) {
  red <- .qn_lta_reduction(state)
  if (is.null(red)) return(NULL)
  w <- state$weights_vec
  unpack <- function(r) {
    st <- .lta_par_unpack(r[red$src], red$state, red$layout)
    if (red$shared) st$tau_independent <- TRUE
    st
  }
  # Per-row log-likelihood and analytic scores come from one forward-backward
  # pass; the finish asks for both at the same point, so it is cached.
  cache <- NULL
  at <- function(r) {
    if (!is.null(cache) && identical(cache$r, r)) return(cache)
    st <- .lta_par_unpack(r[red$src], red$state, red$layout)
    sc <- .lta_score_matrix(st, X, shared_ok = TRUE)
    cache <<- if (is.null(sc)) list(r = r, ll = rep(NA_real_, nrow(X)), S = NULL)
              else list(r = r, ll = sc$ll,
                        # Shared rows fold onto their one free vector.
                        S = t(rowsum(t(sc$S), red$src)))
    cache
  }
  case_ll <- function(r) at(r)$ll
  scores  <- function(r, idx) at(r)$S[, idx, drop = FALSE]
  prior   <- function(r) .lta_log_prior(unpack(r), X, alpha, marginals)
  if (is.na(prior(red$r0))) return(NULL)
  start <- sum(w * case_ll(red$r0)) + prior(red$r0)
  if (!is.finite(start)) return(NULL)
  res <- .qn_finish(case_ll, w, red$r0, prior, scores)
  if (is.null(res) || res$value <= start) return(NULL)

  # New parameters need their own E-step before anything reads a posterior off
  # them; `max_iter = 0` runs exactly that and no EM iteration.
  out <- .lta_em(unpack(res$par), X, max_iter = 0L, alpha = alpha)
  out$converged <- isTRUE(state$converged) || res$converged
  out$n_iter    <- state$n_iter
  out$refined_lbfgs <- state$refined_lbfgs
  out$qn_finish <- list(iterations = res$iterations, n_held = res$n_held,
                        max_gradient = res$max_gradient, converged = res$converged,
                        gain = res$value - start)
  out
}
