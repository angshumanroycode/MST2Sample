\name{euclid_dist}
\alias{euclid_dist}
\alias{mst_edges}
\alias{perm_labels}
\alias{LI_stats}
\alias{I_moments}
\alias{attained_pvalues}
\title{Computational kernels of the pruned-MST tests (internal)}
\description{
\code{euclid_dist} gives the Euclidean distance matrix of the rows of
\code{Z} (as \code{as.matrix(dist(Z))}); \code{mst_edges} the
\eqn{(N-1) \times 2} edge matrix of the MST of a distance matrix (Prim's
algorithm); \code{perm_labels} the matrix of \code{B} random relabelings that
\code{replicate(B, sample(c(rep(1, n), rep(2, m))))} gives on the same seed;
\code{LI_stats} the statistics \eqn{(L_1, L_2, I_1, I_2)} for each
column of a matrix of labels 1/2; \code{I_moments} the exact permutation mean
and covariance matrix of \eqn{(I_1, I_2)} on a graph; and
\code{attained_pvalues} the exact probabilities \eqn{\pi^S} and \eqn{\pi^T}
that a uniformly random relabeling gives \eqn{S = 1} and \eqn{T = 0} on a tree.
}
\usage{
euclid_dist(Z)
mst_edges(D)
perm_labels(n, m, B)
LI_stats(edges, labs)
I_moments(edges, n, m)
attained_pvalues(edges, n, m)
}
\arguments{
\item{Z}{numeric matrix, one observation per row.}
\item{D}{distance matrix.}
\item{edges}{edge matrix with 1-based vertex indices.}
\item{labs}{\eqn{N \times B} matrix of labels 1 and 2, one labeling per column.}
\item{n, m}{sample sizes.}
\item{B}{number of relabelings.}
}
\value{
\code{euclid_dist}: a matrix; \code{mst_edges}: an integer matrix;
\code{LI_stats}: a \eqn{4 \times B} matrix with rows \code{L1}, \code{L2},
\code{I1}, \code{I2}; \code{I_moments}: a list with \code{theta} and \code{V};
\code{attained_pvalues}: a vector with \code{pi_S} and \code{pi_T}.
}
\keyword{internal}
