# Checks of the C++ kernels against direct R computations on small trees.
library(MST2sample)
# the kernels are internal (not exported) from version 0.3.0
for (f in c("euclid_dist", "mst_edges", "perm_labels", "LI_stats", "I_moments",
            "attained_pvalues")) assign(f, getFromNamespace(f, "MST2sample"))

# random tree on N vertices (random attachment)
rtree <- function(N) cbind(sapply(2:N, function(v) sample(v - 1, 1)), 2:N)

components_R <- function(N, E) {          # component id of every vertex
  comp <- seq_len(N)
  repeat {
    changed <- FALSE
    for (e in seq_len(nrow(E))) {
      a <- comp[E[e, 1]]; c <- comp[E[e, 2]]
      if (a != c) { comp[comp == max(a, c)] <- min(a, c); changed <- TRUE }
    }
    if (!changed) return(comp)
  }
}
LI_R <- function(lab, E) {
  N <- length(lab)
  keep <- lab[E[, 1]] == lab[E[, 2]]
  Ek <- E[keep, , drop = FALSE]
  comp <- components_R(N, Ek)
  iso <- tabulate(c(Ek), N) == 0
  c(max(tabulate(comp[lab == 1])), max(tabulate(comp[lab == 2])),
    sum(iso & lab == 1), sum(iso & lab == 2))
}

set.seed(2)
for (rep in 1:30) {
  n <- sample(1:6, 1); m <- sample(1:6, 1); N <- n + m
  E <- rtree(N)
  # all labelings with n ones
  sets <- combn(N, n)
  labs <- apply(sets, 2, function(s) { l <- rep(2, N); l[s] <- 1; l })
  LI <- LI_stats(E, labs)
  ref <- apply(labs, 2, LI_R, E = E)
  stopifnot(all(LI == ref))
  # exact moments of I by enumeration
  mI <- I_moments(E, n, m)
  I <- LI[3:4, , drop = FALSE]
  stopifnot(isTRUE(all.equal(mI$theta, rowMeans(I), check.attributes = FALSE)))
  Vref <- tcrossprod(I - rowMeans(I)) / ncol(I)
  stopifnot(isTRUE(all.equal(mI$V, Vref, check.attributes = FALSE)))
  # exact attained p-values by enumeration
  S <- pmax(LI[1, ] / n, LI[2, ] / m); T <- pmin(LI[3, ] / n, LI[4, ] / m)
  pi <- attained_pvalues(E, n, m)
  stopifnot(isTRUE(all.equal(unname(pi), c(mean(S == 1), mean(T == 0)))))
}

# distance matrix and MST
Z <- matrix(rnorm(15 * 7), 15)
D <- euclid_dist(Z)
stopifnot(identical(D, unname(as.matrix(dist(Z)))))
E <- mst_edges(D)
stopifnot(nrow(E) == 14, all(components_R(15, E) == 1))
# MST weight is minimal among 2000 random spanning trees is not a proof; check
# the cut property instead: every non-tree edge is no shorter than the
# longest edge on the tree path between its ends
adj <- matrix(FALSE, 15, 15); adj[E] <- TRUE; adj[E[, 2:1]] <- TRUE
path_max <- function(a, c) {             # longest edge on the tree path a..c
  prev <- rep(NA, 15); prev[a] <- 0; q <- a
  while (length(q)) { u <- q[1]; q <- q[-1]
    for (v in which(adj[u, ])) if (is.na(prev[v])) { prev[v] <- u; q <- c(q, v) } }
  w <- 0; v <- c
  while (v != a) { w <- max(w, D[v, prev[v]]); v <- prev[v] }
  w
}
for (a in 1:14) for (c in (a + 1):15) if (!adj[a, c])
  stopifnot(D[a, c] >= path_max(a, c))

# path and star (Corollary S.3)
Ep <- cbind(1:39, 2:40)
stopifnot(isTRUE(all.equal(attained_pvalues(Ep, 20, 20)[["pi_S"]], 40 / choose(40, 20))))
Es <- cbind(1, 2:40)
stopifnot(attained_pvalues(Es, 20, 20)[["pi_S"]] == 1)
cat("all kernel checks passed\n")

# perm_labels uses the random numbers exactly as replicate(B, sample(lab))
for (nm in list(c(20, 20), c(20, 40), c(3, 1))) {
  lab <- rep(1:2, nm)
  set.seed(11); a <- replicate(50, sample(lab)); ra <- runif(1)
  set.seed(11); p <- perm_labels(nm[1], nm[2], 50); rp <- runif(1)
  stopifnot(all(a == p), ra == rp)
}
cat("perm_labels matches sample()\n")

# mst_tree: with distinct distances it is mst_edges and uses no random numbers
set.seed(5); Dd <- euclid_dist(matrix(rnorm(30 * 4), 30))
set.seed(6); Et <- mst_tree(Dd); r1 <- runif(1)
set.seed(6); r2 <- runif(1)
stopifnot(identical(Et, mst_edges(Dd)), r1 == r2)

# tied distances (manuscript, Remark "Ties among the pairwise distances"):
# all N = 100 points equal, n = 2, m = 98, b = 99. Without the reordering the
# MST is the star centred at X_1 and T_m rejects at level 0.05 with
# probability about 0.951; with it the level is at most 0.05.
Dz <- matrix(0, 100, 100)
stopifnot(all(mst_edges(Dz)[, 1] == 1))
set.seed(7)
rej <- replicate(400, mst_test_dist(Dz, 2, b = 99, b1 = 2,
                                    ties = "conservative")$p.value[["Tm"]] <= 0.05)
cat("tied distances: rejection rate of T_m =", mean(rej), "\n")
stopifnot(mean(rej) <= 0.1)
