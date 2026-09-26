# Items with different numbers of categories in one "categorical" block. The
# block stores every item at its widest item's width, and until each item
# carried its own count, a two-category item next to a four-category one was
# charged for four: the likelihood was right and n_params, AIC, BIC and every
# likelihood-ratio df were not.

.ragged_fixture <- function(n = 500, seed = 1) {
  set.seed(seed)
  z <- sample(1:2, n, TRUE)
  pick <- function(p1, p2)
    vapply(z, function(k) sample(seq_along(p1), 1L, prob = if (k == 1) p1 else p2),
           integer(1))
  data.frame(
    b1 = pick(c(.8, .2), c(.2, .8)),
    t1 = pick(c(.6, .3, .1), c(.1, .3, .6)),
    q1 = pick(c(.5, .3, .1, .1), c(.1, .1, .3, .5)),
    q2 = pick(c(.4, .4, .1, .1), c(.1, .1, .4, .4)))
}

.fit_cat <- function(d, measurement = "categorical", ...)
  suppressMessages(fit_mixture(d, n_classes = 2, measurement = measurement,
                               n_init = 3, random_state = 1, ...))

test_that("each item is charged for its own categories", {
  d   <- .ragged_fixture()
  fit <- .fit_cat(d)
  expect_equal(fit$mm$cats, c(2L, 3L, 4L, 4L))
  # K * sum(C_j - 1) response probabilities plus K - 1 class weights.
  expect_equal(fit$metrics$n_params, 2 * (1 + 2 + 3 + 3) + 1)
})

test_that("one block, a block per level count and a binary block are one model", {
  d   <- .ragged_fixture()
  one <- .fit_cat(d)
  spl <- .fit_cat(d, list(categorical = "b1", categorical = "t1",
                          categorical = c("q1", "q2")))
  bin <- .fit_cat(d, list(binary = "b1", categorical = c("t1", "q1", "q2")))
  expect_equal(one$metrics$ll, spl$metrics$ll, tolerance = 1e-8)
  expect_equal(one$metrics$ll, bin$metrics$ll, tolerance = 1e-8)
  expect_equal(one$metrics$n_params, spl$metrics$n_params)
  expect_equal(one$metrics$n_params, bin$metrics$n_params)
})

test_that("the summary shows an item's own categories and nothing else", {
  fit <- .fit_cat(.ragged_fixture())
  out <- utils::capture.output(ms <- measurement_summary(fit))
  per_item <- table(ms$item) / 2
  expect_equal(as.vector(per_item[c("b1", "t1", "q1", "q2")]), c(2, 3, 4, 4))
  expect_false(any(grepl("b1 \\(Cat 3\\)", out)))
})

test_that("an unused factor level stays part of the item, with a warning", {
  d <- .ragged_fixture()
  d$q1 <- factor(d$q1, levels = 1:5)
  expect_warning(fit <- .fit_cat(d), "q1 \\(category 5\\)")
  expect_equal(fit$mm$cats, c(2L, 3L, 5L, 4L))
  expect_equal(fit$metrics$n_params, 2 * (1 + 2 + 4 + 3) + 1)
})

test_that("a category empty in one group is kept and reported only where it is free", {
  d <- .ragged_fixture(n = 600, seed = 2)
  g <- factor(rep(c("A", "B"), 300))
  d$q1[g == "B" & d$q1 == 4] <- 3

  expect_warning(
    conf <- .fit_cat(d, group = g, group_effects = "measurement"),
    "q1 in group B \\(category 4\\)")
  # The response space is the pooled one in both groups.
  for (m in conf$mm$models) expect_equal(m$cats, c(2L, 3L, 4L, 4L))
  expect_equal(conf$metrics$n_params, 2 * 2 * (1 + 2 + 3 + 3) + 1)

  # Held equal across groups, the item pools both groups' responses.
  expect_no_warning(
    inv <- .fit_cat(d, group = g, group_effects = "measurement",
                    group_invariant_items = "q1"))
  expect_equal(inv$metrics$n_params, 2 * (1 + 2 + 3) * 2 + 2 * 3 + 1)
})

test_that("a code above an item's own count is refused, not scored as -Inf", {
  X <- cbind(c(1, 2, 3), c(1, 2, 2))
  expect_error(one_hot(X, 3, cats = c(3L, 2L)), NA)
  expect_error(one_hot(cbind(c(1, 2, 3), c(1, 3, 2)), 3, cats = c(3L, 2L)),
               "item 2 has 2 categories")
})
