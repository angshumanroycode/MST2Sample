\name{mst_test}
\alias{mst_test}
\alias{mst_test_dist}
\alias{mst_tree}
\alias{perm_pvalue}
\title{Pruned-MST two-sample tests}
\description{
The tests \eqn{S}, \eqn{T}, \eqn{S_m}, \eqn{T_m} and \eqn{C_m = S_m + T_m}
based on the Euclidean minimum spanning tree (MST) of the pooled sample with
the cross-sample edges removed. \eqn{L_\ell} is the size of the largest
component and \eqn{I_\ell} the number of isolated vertices of sample
\eqn{\ell}; \eqn{S = \max(L_1/n, L_2/m)}, \eqn{T = \min(I_1/n, I_2/m)}, and
\eqn{S_m}, \eqn{T_m} are the Mahalanobis distances of \eqn{(L_1, L_2)} and
\eqn{(I_1, I_2)} from their permutation means, with Moore-Penrose inverses of
the permutation covariance matrices.
}
\usage{
mst_test(X, Y, b = 200, b1 = 1000, ties = c("random", "conservative"))
mst_test_dist(D, n, b = 200, b1 = 1000, ties = c("random", "conservative"))
perm_pvalue(W, ties = c("random", "conservative"))
mst_tree(D)
}
\arguments{
\item{X, Y}{the two samples, one observation per row.}
\item{D}{an \eqn{N \times N} distance matrix of the pooled sample, the first
  \code{n} rows being sample 1.}
\item{n}{size of sample 1.}
\item{b}{number of random relabelings for the p-values.}
\item{b1}{number of auxiliary relabelings for the mean and covariance of
  \eqn{(L_1, L_2)}; they are drawn before, and independently of, the
  \code{b} reference relabelings.}
\item{ties}{\code{"random"}: \eqn{p = (1 + \#\{(W_k, U_k) > (W_0,
  U_0)\})/(b + 1)} with independent uniforms breaking ties (exact size
  \eqn{\lfloor\alpha(b+1)\rfloor/(b+1)}); \code{"conservative"}:
  \eqn{p = (1 + \#\{W_k \ge W_0\})/(b + 1)}.}
\item{W}{observed value followed by the \code{b} permuted values, oriented
  so that large values are evidence against the null hypothesis.}
}
\details{
\code{mst_tree(D)} returns the edges of the pooled MST. If some pairwise
distances are equal, the MST is not unique, and it is computed after a
uniformly random reordering of the pooled sample, so that the choice among
the MSTs does not depend on the sample labels; this is needed for the
level guarantee. With distinct distances the MST is unique, it is the tree that Prim's algorithm gives on
the original order, and no random numbers are used.
The mean and covariance of \eqn{(I_1, I_2)} are computed exactly (Theorem 4.4
of the manuscript). \eqn{S} rejects for large and \eqn{T} for small values;
the exact attained p-values \eqn{\pi^S} on \eqn{\{S = 1\}} and \eqn{\pi^T} on
\eqn{\{T = 0\}} are computed by the dynamic programs of the supplement.
The random numbers are used in the order of the simulation code: \code{b1}
calls of \code{sample()}, \code{b} calls of \code{sample()}, then
\code{runif(b + 1)} for \eqn{S_m}, \eqn{T_m}, \eqn{C_m}, \eqn{S} and
\eqn{T} in turn.
}
\value{
An object of class \code{"mst_test"}: a list with \code{statistic},
\code{p.value}, \code{exact.p.value} (\eqn{\pi^S} if \eqn{S = 1} and
\eqn{\pi^T} if \eqn{T = 0}, else \code{NA}), \code{pi}, \code{L}, \code{I},
\code{theta_L}, \code{V_L}, \code{theta_I}, \code{V_I} and the MST
\code{edges}.
}
\examples{
set.seed(1)
X <- matrix(rnorm(20 * 500), 20)
Y <- matrix(rnorm(20 * 500, 0.2), 20)

## the five tests from the two samples
mst_test(X, Y)

## the same tests from a pooled distance matrix (here the Manhattan
## distance); the first n rows are sample 1
D <- as.matrix(dist(rbind(X, Y), method = "manhattan"))
res <- mst_test_dist(D, n = 20)
res$statistic
res$p.value

## the pooled MST: N - 1 edges with 1-based vertex indices
E <- mst_tree(D)
head(E)
lab <- rep(1:2, c(20, 20))
sum(lab[E[, 1]] != lab[E[, 2]])   # cross-sample edges, removed by the tests

## permutation p-value of another statistic, calibrated as in the package:
## the Friedman-Rafsky number of cross-sample edges, small values being
## evidence against H0, so W is its negative
cross <- function(l) sum(l[E[, 1]] != l[E[, 2]])
W <- -c(cross(lab), replicate(199, cross(sample(lab))))
perm_pvalue(W)                         # randomized tie-break, exact size
perm_pvalue(W, ties = "conservative")  # (1 + #{W_k >= W_0}) / (b + 1)
}
