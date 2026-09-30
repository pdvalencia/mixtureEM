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
#     the empirical information, costs no more than the gradient, and makes
#     the step along such a direction as small as the data say it should be.
#     The scores come from one E-step wherever the model's families have them
#     (R/qn_scores.R, every fit_mixture() model but the growth models; the LTA
#     score blocks for fit_lta()), and from central differences, at 2p
#     passes over the data, otherwise.
#   * Steps are also capped at one unit in any coordinate, and a coordinate
#     that crosses |15| is taken to be heading for the boundary: it is sent to
#     |25| (a probability of 1e-11, the box R/lta_core.R's polish uses) and held
#     there. At |25| its remaining contribution is below 1e-8 of log-likelihood
#     at any sample size fitted here, and a free coordinate there only makes
#     the curvature singular. Only probabilities and rates have such a
#     boundary; the coordinates of a continuous or count model are tagged by
#     kind (.qn_wall_hit()) and a mean or a variance is never walled.
#
# Measured from each fit's own EM winner, the finish reaches the fully
# converged maximum -- the one long EM runs reach -- on every benchmark fit, in
# seconds rather than the tens of thousands of iterations EM would need
# (RECORDS.md, same entry).
#
# The crawl is not a categorical peculiarity. A Gaussian growth mixture with
# missing waves stopped 0.017 short under the same rule, its growth means off
# in the second decimal, and continuous, count, mixed, outcome and
# multiple-group fits all sat short of the maximum by amounts that moved
# reported parameters in the third or fourth decimal (RECORDS.md, "Part 53").
# So the finish covers every model the package fits: each measurement and
# structural model describes its own parameters (R/qn_pack.R), and the one
# model left to its own optimiser is named in .qn_refused with the reason.

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

# The same central differences, returning the diagonal second differences
# too: the points are the ones the gradient already evaluates, so the
# curvature costs one more evaluation in all.
.qn_fd_grad2 <- function(f, par) {
  h  <- 1e-5 * pmax(1, abs(par))
  f0 <- f(par)
  r  <- vapply(seq_along(par), function(i) {
    a <- par; a[i] <- a[i] + h[i]
    b <- par; b[i] <- b[i] - h[i]
    fa <- f(a); fb <- f(b)
    c((fa - fb) / (2 * h[i]), (fa - 2 * f0 + fb) / h[i]^2)
  }, numeric(2))
  list(g = r[1L, ], h = r[2L, ])
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

# What each coordinate of a packed vector is, which decides whether the wall
# applies to it:
#
#   logit    a log-ratio of probabilities (response probabilities, class
#            sizes, class-membership and categorical-outcome coefficients),
#            whose boundary is at either infinity;
#   log_pos  the log of a rate that can reach zero (a Poisson rate), whose
#            boundary is at minus infinity only;
#   scale    the log of a standard deviation, or of a Cholesky diagonal. Its
#            boundary is a degenerate optimum, not a cell at 0 or 1, and it is
#            never walled: a finish that shrinks one is refused instead (see
#            .qn_scale_collapsed());
#   free     anything else -- a mean, a growth coefficient, a regression slope.
#
# A packing that says nothing is all `logit`, which is what every coordinate of
# the categorical fits the finish was first built for is.
.qn_wall_hit <- function(x, kind) {
  (kind == "logit" & abs(x) > .qn_wall) | (kind == "log_pos" & x < -.qn_wall)
}

# The largest fall in a `scale` coordinate the finish may make: log(10), a
# variance divided by 100. The finish starts at a converged EM solution, and a
# genuine move from there is a small correction; a standard deviation falling
# by an order of magnitude is the likelihood's spike at a collapsed variance,
# which EM had been kept out of. The same 1% ratio .gaussian_boundary()
# (R/gaussian_boundary.R) reads as collapse.
.qn_scale_drop <- log(10)

.qn_scale_collapsed <- function(par0, par, kind) {
  s <- which(kind == "scale")
  length(s) > 0L && any(par[s] < par0[s] - .qn_scale_drop)
}

# The finish. `case_ll(par)` is the per-row log-likelihood, `w` the row
# weights, `prior(par)` the log-prior (0 when there is none), and
# `scores(par, idx)` the unweighted per-row scores for the coordinates `idx`
# (central differences of `case_ll` when not supplied). `kind` tags each
# coordinate (see .qn_wall_hit()); `fixed` names coordinates that are held
# where they are -- the measurement block of a two-step or three-step fit.
# `prior_grad(par)`, when supplied, returns the prior's gradient `g` and
# diagonal second derivatives `h` over every coordinate; central differences
# of `prior` otherwise. Returns NULL when nothing is free or the start is not finite, so every
# caller can fall back to the EM solution.
.qn_finish <- function(case_ll, w, par0, prior = function(p) 0,
                       scores = NULL, gtol_stop = .qn_gtol_stop,
                       maxit = 500L, kind = NULL, fixed = integer(0),
                       prior_grad = NULL) {
  scores <- scores %||% function(p, idx) .qn_fd_scores(case_ll, p, idx)
  kind   <- kind %||% rep("logit", length(par0))
  value  <- function(p) {
    v <- sum(w * case_ll(p)) + prior(p)
    if (is.finite(v)) v else -Inf
  }
  x  <- par0
  fx <- value(x)
  if (!is.finite(fx)) return(NULL)
  sw <- sum(w)
  free <- setdiff(which(is.finite(x) & !.qn_wall_hit(x, kind)), fixed)
  if (!length(free)) return(NULL)

  # Coordinates already past the wall go to the edge, if that does not cost.
  to_edge <- function(x, fx, idx) {
    e <- x
    e[idx] <- sign(x[idx]) * .qn_edge
    fe <- value(e)
    if (fe >= fx) list(x = e, fx = fe) else list(x = x, fx = fx)
  }
  held <- setdiff(seq_along(x), c(free, fixed))
  if (length(held)) { r <- to_edge(x, fx, held); x <- r$x; fx <- r$fx }

  it <- 0L
  g  <- NULL
  best_g <- Inf
  stale  <- 0L
  repeat {
    it <- it + 1L
    S  <- scores(x, free)
    gp <- if (is.null(prior_grad)) {
      fr <- function(v) { p <- x; p[free] <- v; prior(p) }
      .qn_fd_grad2(fr, x[free])
    } else {
      a <- prior_grad(x)
      list(g = a$g[free], h = a$h[free])
    }
    g  <- colSums(S * w) + gp$g
    if (max(abs(g)) < gtol_stop * sw || it > maxit) break
    # Ten steps in a row that gained no more than the rounding and left the
    # gradient no smaller: the point is as good as the arithmetic allows. A
    # measurable gain resets the count, since the max-norm of the gradient is
    # not monotone along the way (a cell reaching the wall reshapes it).
    if (max(abs(g)) < best_g) { best_g <- max(abs(g)); stale <- 0L }
    else if ((stale <- stale + 1L) >= 10L) break

    # The curvature is the data's outer product of scores plus the prior's
    # own, which the scores know nothing about. Most priors here are a few
    # pseudo-observations and add little, but the variance prior's curvature
    # in a log standard deviation grows as s^2 / sigma^2, and a class variance
    # far below the item's marginal -- a ceiling effect -- is where it
    # decides the direction: without it the step points along the data's
    # curvature alone, the line search rejects it, and the finish stalls with
    # the gradient far from zero (a free-variance profile model, RECORDS.md,
    # Part 53). Its diagonal, from the differences above, is enough: every
    # measurement prior is separable across parameters. With the priors off
    # it is exactly zero and the step is unchanged.
    B  <- crossprod(S * sqrt(w))
    diag(B) <- diag(B) + pmax(-gp$h, 0)
    e  <- eigen(B, symmetric = TRUE)
    lam <- pmax(e$values, 1e-10 * max(e$values, 1e-300))
    d  <- as.vector(e$vectors %*% (crossprod(e$vectors, g) / lam))
    big <- max(abs(d))
    if (big > .qn_radius) d <- d * .qn_radius / big

    # A step is judged on its gain summed case by case, not on the difference
    # of two totals. Near the maximum a Newton step along a coordinate with a
    # large curvature -- a Gaussian mean in a big class -- gains ~1e-15, below
    # the rounding of a total of |LL| ~ 1e4 (~2e-12), and differencing the
    # totals then reads it as a loss as often as a gain: the line search halves
    # it to nothing and the iterations stall there with the gradient stuck
    # (RECORDS.md, Part 53: a Gaussian LTA). The per-case differences are small
    # numbers whose sum keeps that gain. Even that sum is exact only to the
    # rounding of its terms, ~eps * sum(w |l_i|), and at the last few Newton
    # steps the true gain is that small: a step within it is taken, since it
    # can cost no more than that, and the gradient is what then says whether
    # the finish is still getting anywhere (below). Without the allowance a
    # three-step fit and the same fit on duplicated data stopped at gradients
    # 1e-8 and 5e-7, their estimates 7e-9 apart.
    lx <- case_ll(x)
    px <- prior(x)
    noise <- 4 * .Machine$double.eps * (sum(w * abs(lx)) + abs(px))
    a <- 1
    repeat {
      xn <- x
      xn[free] <- x[free] + a * d
      ln <- case_ll(xn)
      pn <- prior(xn)
      gain <- sum(w * (ln - lx)) + (pn - px)
      fn <- sum(w * ln) + pn
      if (is.finite(fn) && is.finite(gain) && gain >= -noise) break
      a <- a / 2
      if (a < 1e-10) break
    }
    if (is.finite(gain) && gain > noise) stale <- 0L
    # A step that no longer moves any coordinate beyond its last few digits
    # leaves the point where it is, and so would every step after it.
    if (a < 1e-10 ||
        max(abs(xn - x)) <= 1e-14 * max(1, max(abs(x[free])))) break

    hit <- free[.qn_wall_hit(xn[free], kind[free])]
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
# What the finish does not reach, and why
# ------------------------------------------------------------------------------

# Every model the package fits is either taken to the maximum by the finish or
# named here with the reason it is not. The coverage test reads this list, so
# a class that is neither packable nor listed fails it.
.qn_refused <- c(
  lta_random_intercept = paste(
    "A latent transition model with ordinal items and a freely loading",
    "random intercept keeps its own staged search and L-BFGS polish. With",
    "binary items the finish lands on the same optimum and is faster, but",
    "with ordinal items it was slower on every benchmark and landed below",
    "the optimum on one: that optimum has thresholds beyond the logit wall,",
    "which the finish holds at the edge (RECORDS.md, Part 55)."))

# ------------------------------------------------------------------------------
# fit_mixture()
# ------------------------------------------------------------------------------

# The packed vector of a fit and how to read it back: the measurement model,
# then the pooled class sizes unless the structural model supplies the class
# probabilities itself, then the structural model. NULL where any part has no
# packing (R/qn_pack.R).
.qn_mixture_layout <- function(model_state, Y) {
  has_sm   <- !is.null(Y) && !is.null(model_state$sm)
  sm_probs <- has_sm && .supplies_class_probs(model_state$sm)
  K  <- model_state$n_components
  pm <- .qn_pack_mm(model_state$mm)
  if (is.null(pm)) return(NULL)
  parts <- list(mm = pm)
  if (!sm_probs) {
    wts <- pmax(model_state$weights, 1e-12)
    parts$weights <- .qn_part(log(wts[-K] / wts[K]), "logit")
  }
  if (has_sm) {
    ps <- .qn_pack_sm(model_state$sm)
    if (is.null(ps)) return(NULL)
    parts$sm <- ps
  }
  len <- vapply(parts, .qn_len, integer(1))
  end <- cumsum(len)
  idx <- Map(function(e, l) seq_len(l) + e - l, end, len)
  pk  <- .qn_join(parts)
  list(par = pk$par, kind = pk$kind, idx = idx, has_sm = has_sm,
       sm_probs = sm_probs)
}

.qn_mixture_unpack <- function(model_state, lay, par) {
  st <- model_state
  # A frozen block is carried over as it is, not rebuilt from its packing,
  # whose round trip is exact only to the last bit.
  frozen <- model_state$frozen
  if (!"mm" %in% frozen) st$mm <- .qn_unpack_mm(st$mm, par[lay$idx$mm])
  if (!is.null(lay$idx$weights) && !"weights" %in% frozen) {
    lr <- c(par[lay$idx$weights], 0)
    lr <- lr - max(lr)
    st$weights <- exp(lr) / sum(exp(lr))
  }
  if (lay$has_sm) st$sm <- .qn_unpack_sm(st$sm, par[lay$idx$sm])
  st
}

.qn_mixture_supported <- function(model_state, Y) {
  model_state$n_components >= 2L &&
    !is.null(.qn_mixture_layout(model_state, Y))
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

# The objective the finish climbs, as functions of the packed vector: the
# per-case log-likelihood and the log-prior EM's own M-steps add.
.qn_mixture_problem <- function(model_state, X, Y) {
  lay <- .qn_mixture_layout(model_state, Y)
  if (is.null(lay)) return(NULL)
  Yp <- if (lay$has_sm) Y else NULL
  w  <- model_state$sample_weights
  wt <- if (!is.null(w) && length(w) == nrow(X) && any(w != 1)) w else NULL

  # A frozen block (the two-step estimator's measurement model and pooled
  # sizes, R/stepwise.R) stays exactly where it is.
  frozen <- model_state$frozen
  fixed  <- c(if ("mm" %in% frozen) lay$idx$mm,
              if ("weights" %in% frozen) lay$idx$weights)

  unpack <- function(par) .qn_mixture_unpack(model_state, lay, par)
  marg   <- tryCatch(.em_prior_marginals(model_state, X), error = function(e) NULL)
  smp    <- if (lay$has_sm) .qn_sm_prior(model_state$sm, Y, wt) else
    function(sm) 0
  case_ll <- function(par)
    logsumexp(.qn_mixture_logdens(unpack(par), lay, X, Yp), MARGIN = 1)
  prior <- function(par) {
    st <- unpack(par)
    lp <- .em_log_prior(st, X, Yp, marg)
    if (is.na(lp)) NA_real_ else lp + smp(st$sm)
  }
  list(lay = lay, par = lay$par, kind = lay$kind, fixed = fixed, w = w,
       wt = wt, Yp = Yp, unpack = unpack, case_ll = case_ll, prior = prior,
       value = function(par) sum(w * case_ll(par)) + prior(par))
}

.qn_finish_mixture <- function(model_state, X, Y) {
  pr <- .qn_mixture_problem(model_state, X, Y)
  if (is.null(pr)) return(NULL)
  lay <- pr$lay; par0 <- pr$par; Yp <- pr$Yp; wt <- pr$wt
  if (is.na(pr$prior(par0))) return(NULL)
  start <- pr$value(par0)
  if (!is.finite(start)) return(NULL)
  # Analytic scores (R/qn_scores.R) wherever every block has them, from the
  # same E-step as the case log-likelihood: the finish asks for both at one
  # point, so that pass is cached, and the scores are built only when asked
  # for, not at every trial point of the line search. Central differences for
  # the rest -- the growth models.
  s0 <- .qn_mixture_scores(pr$unpack(par0), lay, X, Yp)
  if (!is.null(s0)) {
    cache <- NULL
    at <- function(p) {
      if (is.null(cache) || !identical(cache$p, p)) {
        st <- pr$unpack(p)
        e  <- .qn_mixture_estep(st, lay, X, Yp)
        cache <<- list(p = p, st = st, ll = e$ll, R = e$R, S = NULL)
      }
      cache
    }
    case_ll <- function(p) at(p)$ll
    scores  <- function(p, idx) {
      if (is.null(at(p)$S))
        cache$S <<- .qn_mixture_score_matrix(cache$st, lay, X, Yp, cache$R)
      cache$S[, idx, drop = FALSE]
    }
  } else {
    case_ll <- pr$case_ll
    scores  <- NULL
  }
  # The prior's gradient likewise, where every block has one: differencing it
  # re-derives the data marginals 2p + 1 times an iteration.
  pg <- .qn_mixture_prior_grad(model_state, lay, X, Y, wt)
  if (!is.null(pg) && is.null(pg(pr$unpack(par0)))) pg <- NULL
  prior_grad <- if (!is.null(pg)) function(p) pg(pr$unpack(p))
  # The stopping rule stays the differenced route's, which is the rule every
  # mixture target was graded at (RECORDS.md, "Part 53 graded"): the tighter
  # analytic rule took m_free (group prevalence, p = 76) from 33 iterations to
  # 340 for 2e-8 of log-likelihood (RECORDS.md, Part 56).
  res <- .qn_finish(case_ll, pr$w, par0, pr$prior, scores,
                    gtol_stop = .qn_gtol_stop_fd, kind = lay$kind,
                    fixed = pr$fixed, prior_grad = prior_grad)
  if (is.null(res) || res$value <= start) return(NULL)
  if (.qn_scale_collapsed(par0, res$par, lay$kind)) return(NULL)

  st <- pr$unpack(res$par)
  e  <- e_step(st, X, Yp)
  # An outcome model keeps quantities its M-step derives from the estimates
  # (standard errors, the information matrix, the Wald covariance). One M-step
  # at the finished posteriors rebuilds them there: the point is stationary,
  # so the estimates themselves stay put to within the M-step's own tolerance.
  if (lay$has_sm && !lay$sm_probs) {
    resp  <- exp(e$log_resp)
    st$sm <- if (is.null(wt)) m_step(st$sm, Y, resp) else
      m_step(st$sm, Y, resp, weights = wt)
    e <- e_step(st, X, Yp)
  }
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

# Every latent transition model except one with ordinal items and a freely
# loading random intercept (.qn_refused). Two routes to the scores:
#
#   * analytic, from .lta_score_matrix(), wherever its blocks describe the
#     whole packed vector: binary, ordinal and Gaussian items, and a fixed
#     emission (the three-step's step 3), whose vector is the structural
#     block alone -- one forward-backward pass per iteration;
#   * central differences of .lta_ll_case() for the rest: nominal and count
#     items, a mixed item block, and a step-3 outcome the score blocks do not
#     cover.
#
# A measurement model held at its step-1 estimate (the two-step's step 2) is
# packed with the rest and its coordinates are `fixed`; the model itself is
# carried over untouched, since the packing's round trip is exact only to the
# last bit.
.qn_lta_analytic_families <- c("bernoulli", "bernoulli_nan", "ordinal",
                               "ordinal_nan", "gaussian_diag",
                               "gaussian_diag_nan", "gaussian_unit",
                               "gaussian_unit_nan")

.qn_lta_measurement_fixed <- function(state)
  isTRUE(state$mm_fixed) || "mm" %in% state$frozen

.qn_lta_analytic <- function(state) {
  .lta_scores_supported(state, shared_ok = TRUE) &&
    (isTRUE(state$mm_fixed) ||
       class(state$mm$models[[1]])[1] %in% .qn_lta_analytic_families)
}

# A freely loading random intercept is finished on binary items only. Its
# survivors then run EM to the ordinary rule rather than 1e-11, and skip the
# L-BFGS polish (fit_lta()): the finish lands where the long EM run lands,
# from the same states, in a fraction of the time, measured on seven
# benchmark fits. Ordinal items are refused (.qn_refused). The analytic
# scores cover every coordinate of a binary-item random intercept -- the
# loadings, the node masses, the factor regression -- and nothing else may
# stand in for them: the finite-difference route packs the measurement model
# without the factor.
.qn_lta_ri_finished <- function(state) {
  !is.null(state$ri) && is.null(state$ri$theta) &&
    class(state$mm$models[[1]])[1] %in% c("bernoulli", "bernoulli_nan") &&
    .qn_lta_analytic(state)
}

.qn_lta_supported <- function(state) {
  if (!is.null(state$ri) && .lta_ri_loading_free(state))
    return(.qn_lta_ri_finished(state))
  if (.qn_lta_analytic(state)) return(TRUE)
  if (!is.null(state$ri)) return(FALSE)
  .qn_lta_measurement_fixed(state) || !is.null(.qn_pack_mm(state$mm))
}

# What each block of .lta_par_layout()'s vector is (see .qn_wall_hit()).
# Everything the finish covered before Gaussian items reached it stays a
# logit, so those fits finish exactly as they did. A random intercept's
# loadings, its item-bias shifts and its regression on covariates are
# coefficients with no boundary to head for, so they are never walled.
.qn_lta_block_kind <- function(kind)
  switch(kind, mu = , distal_mu = , lambda = , dif = , ri_beta = "free",
         log_sd = , distal_log_sd = "scale", "logit")

# The packed vector of an LTA state, reduced to the parameters it actually
# has. The transition-free model (`tau_independent`) ties every origin row of
# an occasion's matrix to one vector; with a plain `tau` the ordinary packing
# describes those rows one by one, so they are packed with the flag off and
# mapped back onto one free vector per occasion -- the reduction
# .lta_threestep_vcov() uses for the same model. With a transition regression
# the origin-free design already has one row per occasion and packs as it is.
# `structural_only` packs the structural blocks (and a step-3 distal outcome)
# alone, which is .lta_par_layout()'s own reading of a fixed emission.
.qn_lta_reduction <- function(state, structural_only = FALSE) {
  shared <- isTRUE(state$tau_independent) && is.null(state$tau_beta)
  st <- state
  if (shared) st$tau_independent <- FALSE
  lay_state <- st
  if (structural_only) lay_state$mm_fixed <- TRUE
  lay  <- .lta_par_layout(lay_state)
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
  if (max(abs(r0[src] - full), 0) > 1e-8) return(NULL)
  kind_full <- rep(vapply(lay, function(b) .qn_lta_block_kind(b$kind),
                          character(1)), len)
  list(state = st, layout = lay, src = src, r0 = r0, shared = shared,
       kind = kind_full[match(seq_len(n_red), src)])
}

# The objective's prior. A fixed emission's measurement term is a constant, so
# a step-3 fit carries the structural term alone; everything else is
# .lta_log_prior()'s.
.qn_lta_prior <- function(state, X, alpha, marginals) {
  if (isTRUE(state$mm_fixed)) return(.lta_log_prior_dt(state, alpha))
  .lta_log_prior(state, X, alpha, marginals)
}

.qn_finish_lta <- function(state, X, alpha, marginals) {
  frozen <- "mm" %in% state$frozen
  red <- if (.qn_lta_analytic(state)) .qn_lta_reduction(state) else NULL
  # The analytic route needs the score blocks to describe every coordinate.
  if (!is.null(red)) {
    sc0 <- .lta_score_matrix(state, X, shared_ok = TRUE)
    if (is.null(sc0) || ncol(sc0$S) != length(red$src)) red <- NULL
  }
  analytic <- !is.null(red)
  if (!analytic) red <- .qn_lta_reduction(state, structural_only = TRUE)
  if (is.null(red)) return(NULL)
  mm_free <- !analytic && !.qn_lta_measurement_fixed(state)
  pm <- if (mm_free) .qn_pack_mm(state$mm) else NULL
  if (mm_free && is.null(pm)) return(NULL)
  n_red <- length(red$r0)
  r0    <- c(red$r0, pm$par)
  kind  <- c(red$kind, pm$kind)
  fixed <- if (analytic && frozen)
    unique(red$src[.lta_par_split(red$layout)$measurement]) else integer(0)
  w <- state$weights_vec

  keep <- function(st) {
    if (frozen) {
      st$mm$models <- state$mm$models
      if (!is.null(state$ri)) st$ri <- state$ri
    }
    st
  }
  with_mm <- function(r) {
    st <- red$state
    if (mm_free) st$mm <- .qn_unpack_mm(st$mm, r[-seq_len(n_red)])
    st
  }
  unpack <- function(r) {
    st <- .lta_par_unpack(r[seq_len(n_red)][red$src], with_mm(r), red$layout)
    if (red$shared) st$tau_independent <- TRUE
    keep(st)
  }
  if (analytic) {
    # Per-row log-likelihood and analytic scores come from one forward-backward
    # pass; the finish asks for both at the same point, so it is cached.
    cache <- NULL
    at <- function(r) {
      if (!is.null(cache) && identical(cache$r, r)) return(cache)
      st <- keep(.lta_par_unpack(r[red$src], red$state, red$layout))
      sc <- .lta_score_matrix(st, X, shared_ok = TRUE)
      cache <<- if (is.null(sc)) list(r = r, ll = rep(NA_real_, nrow(X)), S = NULL)
                else list(r = r, ll = sc$ll,
                          # Shared rows fold onto their one free vector.
                          S = t(rowsum(t(sc$S), red$src)))
      cache
    }
    case_ll <- function(r) at(r)$ll
    scores  <- function(r, idx) at(r)$S[, idx, drop = FALSE]
    gtol    <- .qn_gtol_stop
  } else {
    case_ll <- function(r)
      .lta_ll_case(with_mm(r), X, r[seq_len(n_red)][red$src], red$layout)
    scores  <- NULL
    gtol    <- .qn_gtol_stop_fd
  }
  prior <- function(r) .qn_lta_prior(unpack(r), X, alpha, marginals)
  if (is.na(prior(r0))) return(NULL)
  start <- sum(w * case_ll(r0)) + prior(r0)
  if (!is.finite(start)) return(NULL)
  res <- .qn_finish(case_ll, w, r0, prior, scores, gtol_stop = gtol, kind = kind,
                    fixed = fixed)
  if (is.null(res) || res$value <= start) return(NULL)
  if (.qn_scale_collapsed(r0, res$par, kind)) return(NULL)

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
