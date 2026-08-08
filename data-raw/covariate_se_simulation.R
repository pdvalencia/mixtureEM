# ==============================================================================
# Coverage study for the step-three covariate standard errors
# ==============================================================================
#
# Not a data-preparation script (this directory is build-ignored, so it is
# simply where reproducible off-line studies live). It answers the one question
# the reference comparisons in the internal validation suite cannot: a standard
# error can agree with another implementation's and still be the wrong number.
# Only repeated sampling says whether an interval covers.
#
# The design is that of Bakk, Oberski & Vermunt (2014, sec. 4.1), so the results
# are comparable with their Tables 3 and 4: three classes measured by six
# dichotomous indicators, regressed on three covariates each uniform on 1..5,
# class 1 as reference, with
#
#   beta_1. = (-2, 1)   effect of Z1 on classes 2 and 3
#   beta_2. = ( 1, 0)   effect of Z2
#   beta_3. = ( 0, 0)   effect of Z3
#
# and intercepts chosen to give equal class sizes at the covariate means. `rho`,
# the probability of the modal response within a class, sets the separation:
# .80 gives entropy R^2 near .65, .90 near .90.
#
# Four estimators are compared. `q_hessian` is what the package reported before
# R/step3_variance.R existed - the M-step's Q-function Hessian; the other three
# are the values of `se` in fit_mixture().
#
# Run from the package root, as the other scripts here are:
#   Rscript data-raw/covariate_se_simulation.R 250
# Takes roughly an hour at 250 replications; each replication is three full
# three-step fits.

suppressMessages(devtools::load_all("."))

BETA <- rbind(c(0, 0, 0, 0),          # class 1: reference
              c(3, -2, 1, 0),         # intercept, Z1, Z2, Z3
              c(-3, 1, 0, 0))
TRUE_PARS <- c(b12 = -2, b13 = 1, b22 = 1, b23 = 0, b32 = 0, b33 = 0)

gen <- function(n, rho) {
  Z   <- matrix(sample(1:5, n * 3, TRUE), n, 3)
  P   <- softmax_rows(cbind(1, Z) %*% t(BETA))
  cls <- apply(P, 1, function(p) sample(3, 1, prob = p))
  # class 1 positive on all six items; class 2 on the first three; class 3 none
  hi <- rbind(rep(TRUE, 6), c(rep(TRUE, 3), rep(FALSE, 3)), rep(FALSE, 6))
  list(X = matrix(rbinom(n * 6, 1, ifelse(hi[cls, ], rho, 1 - rho)), n, 6),
       Z = Z)
}

# Label switching: match the fitted classes to the generating ones by their
# item-probability profile, since order_by_size is off and class 3 is not
# always the smallest.
align <- function(pis, rho) {
  target <- rbind(rep(rho, 6), c(rep(rho, 3), rep(1 - rho, 3)), rep(1 - rho, 6))
  perms  <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  perms[[which.min(vapply(perms,
    function(p) sum((pis[p, ] - target)^2), numeric(1)))]]
}

one_rep <- function(n, rho) {
  d  <- gen(n, rho)
  Zd <- data.frame(Z1 = d$Z[, 1], Z2 = d$Z[, 2], Z3 = d$Z[, 3])
  out <- list()

  for (m in c("hessian", "robust", "corrected")) {
    f <- try(suppressMessages(suppressWarnings(fit_mixture(
      d$X, Y = Zd, n_components = 3, measurement = "binary",
      structural = "covariate", n_steps = 3, correction = "ML",
      n_init = 5, order_by_size = FALSE, se = m))), silent = TRUE)
    if (inherits(f, "try-error")) return(NULL)

    p   <- align(f$mm$parameters$pis, rho)
    B   <- f$sm$parameters$beta[p, , drop = FALSE]
    D   <- ncol(B)
    idx <- as.vector(sapply(p, function(k) ((k - 1) * D + 1):(k * D)))

    grab <- function(S) {
      S <- S[idx, idx, drop = FALSE]
      est <- se <- numeric(6); nm <- character(6); at <- 0
      for (v in 2:4) for (k in 2:3) {        # Z1..Z3 on classes 2 and 3
        at <- at + 1; nm[at] <- sprintf("b%d%d", v - 1, k)
        ic <- (k - 1) * D + v; ir <- v
        est[at] <- B[k, v] - B[1, v]
        se[at]  <- sqrt(max(0, S[ic, ic] + S[ir, ir] - 2 * S[ic, ir]))
      }
      names(est) <- names(se) <- nm
      list(est = est, se = se, entropy = f$step1_metrics$entropy)
    }

    S <- f$sm$parameters$V_robust
    out[[m]] <- grab(if (is.null(S)) pinv(-f$sm$parameters$hessian) else S)
    # The pre-fix variance. Every branch still stores it, so it costs nothing to
    # read off whichever fit is at hand.
    if (m == "corrected")
      out[["q_hessian"]] <- grab(pinv(-f$sm$parameters$hessian))
  }
  out
}

run <- function(n, rho, reps, seed) {
  set.seed(seed)
  res <- Filter(Negate(is.null), lapply(seq_len(reps), function(i) one_rep(n, rho)))
  ord <- names(TRUE_PARS)
  E   <- t(vapply(res, function(r) r$corrected$est[ord], numeric(6)))

  cat(sprintf("\n=== n = %d, rho = %.2f, %d usable reps, mean entropy R2 = %.3f ===\n",
              n, rho, length(res),
              mean(vapply(res, function(r) r$corrected$entropy, 0))))
  cat("bias:", paste(sprintf("%s %+.3f", ord, colMeans(E) - TRUE_PARS),
                     collapse = "  "), "\n\n")
  cat(sprintf("%-10s %-8s %7s %7s %7s %9s\n",
              "estimator", "par", "sd", "mean se", "se/sd", "coverage"))
  for (m in c("q_hessian", "hessian", "robust", "corrected")) {
    S <- t(vapply(res, function(r) r[[m]]$se[ord], numeric(6)))
    for (j in seq_along(ord))
      cat(sprintf("%-10s %-8s %7.3f %7.3f %7.2f %9.3f\n", m, ord[j],
                  sd(E[, j]), mean(S[, j]), mean(S[, j]) / sd(E[, j]),
                  mean(abs(E[, j] - TRUE_PARS[j]) <= 1.96 * S[, j])))
  }
}

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 250
run(500,  0.80, reps, 4001)   # low separation
run(500,  0.90, reps, 4002)   # high separation, moderate n
run(2000, 0.90, reps, 4003)   # high separation, large n
