# Pruned-MST two-sample tests S, T, S_m, T_m and C_m (manuscript_v11.tex,
# Sections 3 and 4, with the calibration of Section 7).

mahal <- function(x, theta, Vinv) {   # columns of x are (L1, L2) or (I1, I2)
  z <- x - theta
  colSums(z * (Vinv %*% z))
}

# Permutation p-value of W[1] against W[-1] (large values are evidence
# against H0). "random": p = (1 + #{(W_k, U_k) > (W_0, U_0)}) / (b + 1) with
# independent uniforms U breaking ties, which has exact size
# floor(alpha (b + 1)) / (b + 1); "conservative": p of equation (4.1),
# (1 + #{W_k >= W_0}) / (b + 1).
perm_pvalue <- function(W, ties = c("random", "conservative")) {
  ties <- match.arg(ties)
  Wk <- W[-1]
  if (ties == "random") {
    U  <- runif(length(W))
    Uk <- U[-1]
    (1 + sum(Wk > W[1] | (Wk == W[1] & Uk > U[1]))) / length(W)
  } else {
    (1 + sum(Wk >= W[1])) / length(W)
  }
}

# Pooled MST. With tied pairwise distances the MST is not unique, and Prim's
# algorithm would choose among the MSTs by the order of the observations,
# i.e. by the sample labels, which breaks the level of the tests
# (manuscript, Remark "Ties among the pairwise distances"). So the MST is then
# computed after a uniformly random reordering of the pooled sample. With
# distinct distances the MST is unique, nothing is reordered and no random
# numbers are used, so earlier results reproduce exactly.
mst_tree <- function(D) {
  D <- as.matrix(D)
  storage.mode(D) <- "double"
  if (!anyDuplicated(D[upper.tri(D)])) return(mst_edges(D))
  o <- sample.int(nrow(D))
  E <- mst_edges(D[o, o, drop = FALSE])
  a <- o[E[, 1]]; c <- o[E[, 2]]
  cbind(pmin(a, c), pmax(a, c))
}

mst_test <- function(X, Y, b = 200, b1 = 1000,
                     ties = c("random", "conservative")) {
  dname <- paste(deparse1(substitute(X)), "and", deparse1(substitute(Y)))
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (ncol(X) != ncol(Y)) stop("X and Y must have the same number of columns")
  storage.mode(X) <- "double"; storage.mode(Y) <- "double"
  D   <- euclid_dist(rbind(X, Y))
  out <- mst_test_dist(D, nrow(X), b = b, b1 = b1, ties = ties)
  out$data.name <- dname
  out
}

mst_test_dist <- function(D, n, b = 200, b1 = 1000,
                          ties = c("random", "conservative")) {
  ties <- match.arg(ties)
  D <- as.matrix(D)
  storage.mode(D) <- "double"
  N <- nrow(D); m <- N - n
  if (ncol(D) != N) stop("D must be a square matrix")
  if (n < 1 || m < 1) stop("both samples must be non-empty")
  if (b1 < 2) stop("b1 must be at least 2")
  lab   <- c(rep(1, n), rep(2, m))
  edges <- mst_tree(D)

  # theta_I, V_I exact (Theorem 4.4); theta_L, V_L from b1 auxiliary
  # relabelings drawn before, and independently of, the b reference ones
  mI    <- I_moments(edges, n, m)
  VIinv <- ginv(mI$V)
  aux   <- LI_stats(edges, perm_labels(n, m, b1))[1:2, , drop = FALSE]
  thL   <- rowMeans(aux)
  VL    <- cov(t(aux))
  VLinv <- ginv(VL)

  # observed labels (column 1) and b reference relabelings, shared by all tests
  labs <- cbind(lab, perm_labels(n, m, b))
  LI   <- LI_stats(edges, labs)
  Sm   <- mahal(LI[1:2, , drop = FALSE], thL, VLinv)
  Tm   <- mahal(LI[3:4, , drop = FALSE], mI$theta, VIinv)
  Cm   <- Sm + Tm
  S    <- pmax(LI[1, ] / n, LI[2, ] / m)
  T    <- pmin(LI[3, ] / n, LI[4, ] / m)

  p  <- c(Sm = perm_pvalue(Sm, ties), Tm = perm_pvalue(Tm, ties),
          Cm = perm_pvalue(Cm, ties), S = perm_pvalue(S, ties),
          T = perm_pvalue(-T, ties))
  pi <- attained_pvalues(edges, n, m)
  structure(list(
    statistic = c(Sm = Sm[1], Tm = Tm[1], Cm = Cm[1], S = S[1], T = T[1]),
    p.value   = p,
    exact.p.value = c(S = if (S[1] == 1) pi[["pi_S"]] else NA_real_,
                      T = if (T[1] == 0) pi[["pi_T"]] else NA_real_),
    pi = pi, L = LI[1:2, 1], I = LI[3:4, 1],
    theta_L = thL, V_L = VL, theta_I = mI$theta, V_I = mI$V,
    edges = edges, n = n, m = m, b = b, b1 = b1, ties = ties,
    data.name = "D"), class = "mst_test")
}

print.mst_test <- function(x, digits = 4, ...) {
  cat("\n\tPruned-MST two-sample tests\n\n")
  cat("data: ", x$data.name, "  (n = ", x$n, ", m = ", x$m, ")\n", sep = "")
  cat("L = (", paste(x$L, collapse = ", "), "), I = (",
      paste(x$I, collapse = ", "), ")\n", sep = "")
  tab <- data.frame(statistic = signif(x$statistic, digits),
                    p.value = signif(x$p.value, digits))
  print(tab)
  cat("\np-values from b = ", x$b, " relabelings (ties broken ",
      if (x$ties == "random") "at random" else "conservatively",
      "); theta_L, V_L from b1 = ", x$b1, " auxiliary relabelings.\n", sep = "")
  cat("Exact attained p-values: pi_S = ", signif(x$pi[["pi_S"]], digits),
      ", pi_T = ", signif(x$pi[["pi_T"]], digits),
      " (the p-values of S on {S = 1} and of T on {T = 0}).\n\n", sep = "")
  invisible(x)
}
