test_that(".gauss_hermite() reduces to a point mass at Q = 1", {
  q <- .gauss_hermite(1)
  expect_equal(q$z, 0)
  expect_equal(q$w, 1)
})

test_that(".gauss_hermite() weights and nodes are correct at Q = 5", {
  q <- .gauss_hermite(5)
  expect_equal(sort(q$z),
               c(-2.8569700139, -1.3556261800, 0, 1.3556261800, 2.8569700139),
               tolerance = 1e-9)
  expect_equal(q$w[order(q$z)],
               c(0.0112574113, 0.2220759220, 0.5333333333, 0.2220759220,
                 0.0112574113),
               tolerance = 1e-9)
})

test_that(".gauss_hermite() reproduces standard-normal moments", {
  for (Q in c(5, 15, 20, 30, 100)) {
    q <- .gauss_hermite(Q)
    expect_equal(sum(q$w), 1, tolerance = 1e-13)
    expect_equal(sum(q$w * q$z), 0, tolerance = 1e-12)
    # This is the check that catches a physicists'-rule Hermite implementation
    # (weight exp(-x^2) rather than the probabilists' exp(-x^2/2)): a
    # physicists' rule gives 0.5 here, not 1, because its nodes want dividing
    # by sqrt(2) before they mean anything on this scale.
    expect_equal(sum(q$w * q$z^2), 1, tolerance = 1e-13)
    expect_equal(sum(q$w * q$z^4), 3, tolerance = 1e-11)
  }
})

test_that(".gauss_hermite() integrates plogis(1 + 2z), and Q = 5 is not enough", {
  # Measured 2026-09-01 (roadmap ### 14.6 item 1). The paper's own 0.644 is a
  # different quantity (its eq.-(4) probit approximation) and is not asserted
  # here -- see internal/ROADMAP.md ### 14.10.1.
  targets <- c(`15` = 0.6477238, `20` = 0.6477294, `30` = 0.6477267,
               `100` = 0.6477264)
  for (Q in c(15, 20, 30, 100)) {
    q <- .gauss_hermite(Q)
    val <- sum(q$w * plogis(1 + 2 * q$z))
    expect_equal(val, targets[[as.character(Q)]], tolerance = 1e-4)
  }

  q5 <- .gauss_hermite(5)
  val5 <- sum(q5$w * plogis(1 + 2 * q5$z))
  expect_equal(val5, 0.6519961, tolerance = 1e-4)
  expect_gt(abs(val5 - 0.6477264), 1e-3)
})
