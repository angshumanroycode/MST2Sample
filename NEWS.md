# MST2sample 0.3.0

* Smaller public interface (2026-11-12). Only `mst_test`, `mst_test_dist`,
  `perm_pvalue` and `mst_tree` are exported. The C++ kernels `euclid_dist`,
  `mst_edges`, `perm_labels`, `LI_stats`, `I_moments` and `attained_pvalues`
  are now internal; call them as `MST2sample:::name`. No computation changes.

* `print.mst_test` is documented on its own internal help page, so the help
  index lists only the four exported functions.

* The help page has examples for `mst_test_dist`, `mst_tree` and `perm_pvalue`.

# MST2sample 0.2.0

* Tie handling in the tree construction (2026-10-24). If some pooled pairwise
  distances tie, `mst_tree` now computes the MST after a uniformly random
  reordering of the pooled sample, so that the choice among the MSTs does not
  depend on the sample labels (needed for the level of the tests). With
  distinct distances nothing changes and no random numbers are used. No
  reported number for the Pomeroy and Gravier examples changes (they have no
  tied distances).

# MST2sample 0.1.0

* First version: S, T, S_m, T_m and C_m, the exact permutation moments of I,
  and the exact attained p-values of S and T, computed in C++.
