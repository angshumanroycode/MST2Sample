# MST2sample

Version 0.3.0 (see `NEWS.md`).

R package (Rcpp) for the pruned-MST two-sample tests S, T, S_m, T_m and C_m of
`manuscript_v11.tex`, calibrated as in Section 7.

## Install

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("angshumanroycode/MST2sample")
```

The package contains C++ code, so a compiler toolchain is needed: Rtools on
Windows, the Xcode command line tools on macOS (`xcode-select --install`), or
`build-essential` (or your distribution's equivalent) on Linux.

Alternatively, install the source tarball `MST2sample_0.3.0.tar.gz`:

```r
install.packages("MST2sample_0.3.0.tar.gz", repos = NULL, type = "source")
```

## Use

```r
library(MST2sample)
res <- mst_test(X, Y, b = 200, b1 = 1000)   # rows are observations
res$statistic; res$p.value; res$pi
mst_test_dist(D, n)                          # from a pooled distance matrix
```

What it computes:

* the Euclidean distance matrix and the pooled MST (Prim, O(N^2)) in C++;
  if some pairwise distances tie, `mst_tree` computes the MST after a uniformly
  random reordering of the pooled sample, so that the choice among the MSTs does
  not depend on the sample labels (needed for the level; with distinct
  distances nothing changes and no random numbers are used);
* L = (L1, L2) and I = (I1, I2) for every relabeling by one union-find pass in C++;
* theta_I and V_I exactly (Theorem 4.4), with the same arithmetic as the R code;
* theta_L and V_L from b1 auxiliary relabelings drawn before the b reference relabelings;
* Moore-Penrose inverses (`MASS::ginv`);
* p = (1 + #{(W_k, U_k) > (W_0, U_0)}) / (b + 1) with a randomized tie-break
  (`ties = "conservative"` gives (1 + #{W_k >= W_0}) / (b + 1));
* S rejects for large and T for small values; the exact attained p-values
  pi_S and pi_T (Proposition S.2) by the dynamic programs of the supplement.

Exported: `mst_test`, `mst_test_dist`, `perm_pvalue` and `mst_tree`. The C++
kernels `euclid_dist`, `mst_edges`, `perm_labels`, `LI_stats`, `I_moments` and
`attained_pvalues` are internal (from 0.3.0); the simulation scripts call them
as `MST2sample:::name`.

## Checks

* `tests/test-kernels.R` (run by `R CMD check`) compares the kernels with full
  enumeration over all labelings on random small trees, and checks pi_S on a
  path and a star (Corollary S.3).
