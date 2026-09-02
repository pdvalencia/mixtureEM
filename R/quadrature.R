# Gauss-Hermite quadrature nodes and weights for a standard normal density,
# via the Golub-Welsch eigendecomposition of the recurrence's Jacobi matrix.
# `Q = 1` is special-cased: its node is 0 and its weight is 1, which is also
# the case that lets a random-intercept fit reduce exactly to a model with no
# random intercept at all (used as a regression test in test-lta-ri.R).
.gauss_hermite <- function(Q) {
  Q <- as.integer(Q)
  if (Q < 1L) stop("`Q` must be at least 1.", call. = FALSE)
  if (Q == 1L) return(list(z = 0, w = 1))
  J  <- matrix(0, Q, Q)
  od <- sqrt(seq_len(Q - 1L))
  for (i in seq_len(Q - 1L)) { J[i, i + 1L] <- od[i]; J[i + 1L, i] <- od[i] }
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(z = e$values[o], w = (e$vectors[1L, o])^2)
}
