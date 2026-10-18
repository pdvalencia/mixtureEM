# predict(): posterior class probabilities for cases outside the estimation
# sample, at the fitted parameters.
#
# Two contracts. On the estimation sample of a one-step model it must equal
# class_assignments() exactly, since both are the same E-step. On held-out
# cases it must reproduce posteriors computed at the same fixed parameters,
# missing answers and a covariate included; those are pinned below to 1e-6.

.p_bin_fit <- function() {
  set.seed(3301)
  P <- rbind(c(.9, .85, .8, .85, .9, .8), c(.2, .25, .8, .85, .3, .2),
             c(.1, .15, .1, .2, .15, .1))
  z <- sample(3, 600, TRUE, c(.45, .35, .20))
  X <- t(sapply(z, function(k) rbinom(6, 1, P[k, ])))
  colnames(X) <- paste0("u", 1:6)
  list(X = X, fit = suppressMessages(
    fit_mixture(X, n_classes = 3, measurement = "binary", n_init = 20,
                random_state = 1, n_cores = 1)))
}

.p_lpa_data <- function(n) {
  x  <- rnorm(n)
  z  <- ifelse(runif(n) < plogis(-0.3 + 1.2 * x), 1L, 2L)
  mu <- rbind(c(0, 0, 0), c(1.8, 1.5, 2.0))
  s  <- rbind(c(1, 1, 1), c(.6, 1.4, .8))
  Y  <- t(sapply(z, function(k) rnorm(3, mu[k, ], s[k, ])))
  d  <- data.frame(Y, x = x)
  names(d)[1:3] <- paste0("y", 1:3)
  d
}

test_that("held-out binary cases, missing answers included, match pinned posteriors", {
  b  <- .p_bin_fit()
  ho <- rbind(c(NA, NA, NA, NA, NA, 1), c(0, 0, 1, 1, 1, 0),
              c(1, 1, 0, 1, 0, 1),      c(1, 0, 0, 0, 1, 0),
              c(0, 1, 0, NA, 0, NA))
  colnames(ho) <- paste0("u", 1:6)
  target <- rbind(c(0.7982852749, 0.0907974072, 0.1109173179),
                  c(0.0161916151, 0.0715498806, 0.9122585043),
                  c(0.9397691658, 0.0140173037, 0.0462135305),
                  c(0.1149103911, 0.8792759481, 0.0058136608),
                  c(0.0284257225, 0.7592819949, 0.2122922826))
  post <- predict(b$fit, ho, type = "posterior")
  expect_lt(max(abs(unname(post) - target)), 1e-6)
  expect_identical(predict(b$fit, ho), max.col(target))
})

test_that("held-out cases of a model with a predictor match pinned posteriors", {
  set.seed(3302)
  d   <- .p_lpa_data(600)
  fit <- suppressMessages(
    fit_mixture(d[, 1:3], n_classes = 2, measurement = "continuous",
                predictors = d["x"], variances_equal = FALSE, n_init = 20,
                random_state = 1, n_cores = 1))
  ho <- data.frame(y1 = c(0.1936675, 2.6329611, 0.5284713, -0.7323558,
                          0.9985243, 0.9517388),
                   y2 = c(1.2023066, 2.4345464, 0.4392447, 1.9413765,
                          -1.0424929, NA),
                   y3 = c(-0.1476771, 2.4989400, -1.4958281, 0.4024714,
                          -1.1508158, NA),
                   x  = c(0.55414317, -0.46107692, 0.08957698, 1.48521138,
                          0.23345658, 0.63520182))
  target <- c(0.0045054073, 0.9999654084, 0.0001053499, 0.0002779340,
              0.0004436977, 0.4549809350)
  # The data above are printed to 7-8 digits, so the pin is looser than the
  # 1e-8 the full-precision values reach.
  expect_lt(max(abs(predict(fit, ho, "posterior")[, 1] - target)), 1e-6)
})

test_that("on the estimation sample predict() equals class_assignments()", {
  b <- .p_bin_fit()
  expect_equal(predict(b$fit, b$X, "posterior"),
               class_assignments(b$fit, "posterior"), tolerance = 1e-12)
  expect_identical(predict(b$fit, b$X), class_assignments(b$fit))

  set.seed(5)
  n  <- 300
  z  <- sample(2, n, TRUE)
  d  <- data.frame(p1 = factor(ifelse(z == 1, sample(c("lo", "mid", "hi"), n, TRUE, c(.7, .2, .1)),
                                      sample(c("lo", "mid", "hi"), n, TRUE, c(.1, .2, .7))),
                               levels = c("lo", "mid", "hi")),
                   p2 = sample(1:3, n, TRUE),
                   c1 = rnorm(n, 2 * z), k1 = rpois(n, 2 * z),
                   b1 = ifelse(runif(n) < .3 * z, "yes", "no"))
  d$c1[c(3, 9)] <- NA
  fits <- suppressWarnings(suppressMessages(list(
    poly = fit_mixture(d[c("p1", "p2")], 2, "categorical", n_init = 3,
                        random_state = 1, n_cores = 1),
    cont  = fit_mixture(d["c1"], 2, "continuous", n_init = 3,
                        random_state = 1, n_cores = 1),
    count = fit_mixture(d["k1"], 2, "count", n_init = 3, random_state = 1,
                        n_cores = 1),
    mixed = fit_mixture(d, 2, list(categorical = 1:2, continuous = 3,
                                   count = 4, binary = 5),
                        n_init = 3, random_state = 1, n_cores = 1),
    cov1  = fit_mixture(d["c1"], 2, "continuous", predictors = d["k1"],
                        n_steps = 1, n_init = 3, random_state = 1,
                        n_cores = 1))))
  # Rows 3 and 9 miss the one continuous item, so the single-item fits dropped
  # them; predict() keeps their rows and gives them NA.
  for (nm in names(fits)) {
    f    <- fits[[nm]]
    keep <- if (nm %in% c("cont", "cov1")) -c(3, 9) else seq_len(n)
    p    <- predict(f, d, "posterior")
    if (nm %in% c("cont", "cov1")) expect_true(all(is.na(p[c(3, 9), ])))
    expect_equal(p[keep, ],
                 class_assignments(f, "posterior"), tolerance = 1e-10,
                 label = nm)
  }
})

test_that("new data are matched by name and coded by the fitted levels", {
  set.seed(6)
  n <- 200
  z <- sample(2, n, TRUE)
  d <- data.frame(a = factor(ifelse(runif(n) < .3 * z, "yes", "no")),
                  b = factor(sample(c("x", "y", "z"), n, TRUE)),
                  c = ifelse(runif(n) < .4 * z, 2, 5))
  fit <- suppressMessages(fit_mixture(d, 2, list(binary = c(1, 3),
                                                 categorical = 2),
                                      n_init = 3, random_state = 1,
                                      n_cores = 1))
  ref <- predict(fit, d, "posterior")
  # Reordered, with an extra column, and with a factor whose levels come in a
  # different order: the same answers must give the same posteriors.
  d2 <- d[, c("c", "b", "a")]
  d2$extra <- 1
  d2$b <- factor(as.character(d2$b), levels = c("z", "y", "x"))
  expect_equal(predict(fit, d2, "posterior"), ref, tolerance = 1e-12)
  # A character column is coded the same way as the factor it replaces.
  d3 <- d
  d3$b <- as.character(d3$b)
  expect_equal(predict(fit, d3, "posterior"), ref, tolerance = 1e-12)

  expect_error(predict(fit, d[, c("a", "b")]), "lacks indicator")
  d4 <- d[1:3, ]
  d4$b <- factor(c("x", "w", "y"))
  expect_error(predict(fit, d4), "never saw: w")
  d5 <- d[1:3, ]
  d5$b <- c(1, 2, 3)
  expect_error(predict(fit, d5), "was a factor")
  d6 <- d[1:3, ]
  d6$c <- c(2, 5, 7)
  expect_error(predict(fit, d6), "other than its two")
})

test_that("unscorable cases get NA and keep their row", {
  b  <- .p_bin_fit()
  ho <- b$X[1:4, ]
  ho[2, ] <- NA
  out <- predict(b$fit, ho, "both")
  expect_identical(nrow(out), 4L)
  expect_true(is.na(out$class[2]) && all(is.na(out[2, -1])))
  expect_false(anyNA(out[-2, ]))

  set.seed(3302)
  d   <- .p_lpa_data(300)
  fit <- suppressMessages(
    fit_mixture(d[, 1:3], 2, "continuous", predictors = d["x"], n_steps = 1,
                n_init = 3, random_state = 1, n_cores = 1))
  nd <- d[1:3, ]
  nd$x[3] <- NA
  expect_warning(p <- predict(fit, nd), "missing predictor")
  expect_true(is.na(p[3]) && !anyNA(p[1:2]))
  expect_error(predict(fit, d[1:3, 1:3]), "lacks predictor 'x'")
})

test_that("a distal outcome never enters a new case's classification", {
  set.seed(8)
  n <- 300
  z <- sample(2, n, TRUE)
  X <- matrix(rbinom(n * 5, 1, ifelse(z == 1, .8, .2)), n,
              dimnames = list(NULL, paste0("u", 1:5)))
  y <- rnorm(n, 3 * z)
  fit <- suppressMessages(fit_mixture(X, 2, "binary", outcome = y,
                                      n_steps = 1, n_init = 3,
                                      random_state = 1, n_cores = 1))
  noy <- fit
  noy$sm <- NULL
  expect_equal(unname(predict(fit, X, "posterior")),
               exp(mixtureEM:::e_step(noy, fit$data, NULL)$log_resp),
               tolerance = 1e-12)
})

test_that("models predict() does not cover are refused", {
  set.seed(22)
  X <- matrix(rbinom(150 * 3, 1, 0.5), ncol = 3)
  lta <- suppressWarnings(
    fit_lta(X, n_statuses = 2, times = 3, measurement = "binary", n_init = 2,
            random_state = 2, standard_errors = FALSE))
  expect_error(predict(lta, X), "latent transition")

  b <- .p_bin_fit()
  expect_error(predict(b$fit), "newdata")
  g <- b$fit
  g$group_info <- list()
  expect_error(predict(g, b$X), "multiple-group")
})
