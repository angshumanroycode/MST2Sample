// Pruned-MST two-sample tests: the computational kernels.
//
// The arithmetic follows the R code in "Codes by Claude/Figure 2 equal
// sample.R" operation by operation (same summation order, long double
// accumulation where R's sum() uses it, Rmath's lgammafn), so that the
// statistics agree with that code to the last bit and the randomized
// permutation p-values agree exactly on the same seed.

#include <Rcpp.h>
#include <Rmath.h>
#include <vector>
#include <limits>
#include <cmath>

using namespace Rcpp;

// ── Euclidean distance matrix, as.matrix(dist(Z)) ───────────────────────────
// [[Rcpp::export]]
NumericMatrix euclid_dist(NumericMatrix Z) {
  const int N = Z.nrow(), d = Z.ncol();
  NumericMatrix D(N, N);
  const double *x = Z.begin();
  for (int i = 0; i < N; ++i)
    for (int j = i + 1; j < N; ++j) {
      double s = 0.0;
      const double *xi = x + i, *xj = x + j;
      for (int k = 0; k < d; ++k, xi += N, xj += N) {
        double dev = *xi - *xj;
        s += dev * dev;
      }
      s = std::sqrt(s);
      D(i, j) = s;
      D(j, i) = s;
    }
  return D;
}

// ── Minimum spanning tree of a dense distance matrix (Prim, O(N^2)) ─────────
// Returns the N - 1 edges as an (N - 1) x 2 matrix of 1-based vertex indices.
// [[Rcpp::export]]
IntegerMatrix mst_edges(NumericMatrix D) {
  const int N = D.nrow();
  if (D.ncol() != N) stop("D must be a square matrix");
  if (N < 2) stop("need at least two points");
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> key(N, inf);
  std::vector<int> from(N, -1);
  std::vector<char> in(N, 0);
  IntegerMatrix E(N - 1, 2);
  int cur = 0;
  in[0] = 1;
  for (int step = 0; step < N - 1; ++step) {
    int best = -1;
    double bk = inf;
    for (int v = 0; v < N; ++v) {
      if (in[v]) continue;
      double w = D(cur, v);
      if (w < key[v]) { key[v] = w; from[v] = cur; }
      if (best < 0 || key[v] < bk) { bk = key[v]; best = v; }
    }
    in[best] = 1;
    int a = from[best], c = best;
    E(step, 0) = (a < c ? a : c) + 1;
    E(step, 1) = (a < c ? c : a) + 1;
    cur = best;
  }
  return E;
}

// ── L = (L1, L2) and I = (I1, I2) for many labelings ────────────────────────
// edges: (N - 1) x 2 (1-based); labs: N x B matrix of labels 1/2, one labeling
// per column. Returns the 4 x B matrix with rows L1, L2, I1, I2. One
// union-find pass over the edges per labeling.
static inline int uf_find(std::vector<int> &p, int x) {
  while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
  return x;
}

// [[Rcpp::export]]
NumericMatrix LI_stats(IntegerMatrix edges, NumericMatrix labs) {
  const int N = labs.nrow(), B = labs.ncol(), nE = edges.nrow();
  std::vector<int> e1(nE), e2(nE);
  for (int e = 0; e < nE; ++e) {
    e1[e] = edges(e, 0) - 1; e2[e] = edges(e, 1) - 1;
    if (e1[e] < 0 || e1[e] >= N || e2[e] < 0 || e2[e] >= N)
      stop("edge index out of range");
  }
  NumericMatrix out(4, B);
  std::vector<int> par(N), sz(N), lab(N);
  std::vector<char> touched(N);
  for (int k = 0; k < B; ++k) {
    for (int v = 0; v < N; ++v) {
      double l = labs(v, k);
      if (l != 1.0 && l != 2.0) stop("labels must be 1 or 2");
      lab[v] = (int) l; par[v] = v; sz[v] = 1; touched[v] = 0;
    }
    for (int e = 0; e < nE; ++e) {
      int a = e1[e], c = e2[e];
      if (lab[a] != lab[c]) continue;          // cross-sample edge removed
      touched[a] = touched[c] = 1;
      int ra = uf_find(par, a), rc = uf_find(par, c);
      if (ra == rc) continue;
      if (sz[ra] < sz[rc]) std::swap(ra, rc);
      par[rc] = ra; sz[ra] += sz[rc];
    }
    int L[3] = {0, 0, 0}, I[3] = {0, 0, 0};
    for (int v = 0; v < N; ++v) {
      if (par[v] == v && sz[v] > L[lab[v]]) L[lab[v]] = sz[v];
      if (!touched[v]) ++I[lab[v]];
    }
    out(0, k) = L[1]; out(1, k) = L[2]; out(2, k) = I[1]; out(3, k) = I[2];
  }
  rownames(out) = CharacterVector::create("L1", "L2", "I1", "I2");
  return out;
}

// ── Random relabelings, as replicate(B, sample(lab)) ────────────────────────
// lab = c(rep(1, n), rep(2, m)). Uses R's random number stream with the
// algorithm of R's sample() without replacement (R >= 3.6.0, "Rejection"
// sampling), so it gives the same matrix as replicate(B, sample(lab)) after
// the same set.seed().
// [[Rcpp::export]]
NumericMatrix perm_labels(int n, int m, int B) {
  const int N = n + m;
  NumericMatrix labs(N, B);
  std::vector<int> x(N);
  RNGScope scope;
  for (int k = 0; k < B; ++k) {
    for (int i = 0; i < N; ++i) x[i] = i;
    int left = N;
    for (int i = 0; i < N; ++i) {
      int j = (int) R_unif_index(left);
      labs(i, k) = x[j] < n ? 1.0 : 2.0;
      x[j] = x[--left];
    }
  }
  return labs;
}

// ── Exact permutation moments of I (Theorem 4.4) ────────────────────────────
// lfall(a, k) = log of the falling factorial (a)_k, and h(k, l) =
// (n)_k (m)_l / (N)_{k+l}, computed exactly as lfall() and h() in the R code.
static inline double lfall(double a, double k) {
  return (k <= a) ? Rf_lgammafn(a + 1) - Rf_lgammafn(a - k + 1)
                  : R_NegInf;
}
static inline double hfun(double k, double l, double n, double m) {
  if (!(k <= n && l <= m)) return 0.0;
  double kk = std::min(k, n), ll = std::min(l, m);
  return std::exp(lfall(n, kk) + lfall(m, ll) - lfall(n + m, kk + ll));
}

// [[Rcpp::export]]
List I_moments(IntegerMatrix edges, int n, int m) {
  const int N = n + m, nE = edges.nrow();
  std::vector<std::vector<int> > nb(N);
  std::vector<char> adj((size_t) N * N, 0);
  for (int e = 0; e < nE; ++e) {
    int a = edges(e, 0) - 1, c = edges(e, 1) - 1;
    if (a < 0 || a >= N || c < 0 || c >= N) stop("edge index out of range");
    if (a == c || adj[(size_t) a * N + c]) continue;   // simple graph
    adj[(size_t) a * N + c] = adj[(size_t) c * N + a] = 1;
    nb[a].push_back(c); nb[c].push_back(a);
  }
  std::vector<double> deg(N);
  for (int v = 0; v < N; ++v) deg[v] = (double) nb[v].size();
  // c_uv = number of common neighbours, from the paths u - w - v
  std::vector<int> C((size_t) N * N, 0);
  for (int w = 0; w < N; ++w)
    for (int u : nb[w]) for (int v : nb[w]) ++C[(size_t) u * N + v];

  const double dn = n, dm = m;
  long double EI1 = 0, EI2 = 0;
  for (int v = 0; v < N; ++v) {
    EI1 += hfun(1, deg[v], dn, dm);
    EI2 += hfun(deg[v], 1, dn, dm);
  }
  // sums over u != v in column-major order (u the row, v the column), as
  // sum(EAA[off]) in R
  long double sAA = 0, sBB = 0, sAB = 0;
  for (int v = 0; v < N; ++v)
    for (int u = 0; u < N; ++u) {
      if (u == v) continue;
      size_t uv = (size_t) u * N + v;
      double du = deg[u], dv = deg[v], c = C[uv];
      if (adj[uv]) {
        sAA += 0.0; sBB += 0.0;
        sAB += hfun(dv, du, dn, dm) * (double) (c == 0);
      } else {
        sAA += hfun(2, du + dv - c, dn, dm);
        sBB += hfun(du + dv - c, 2, dn, dm);
        sAB += hfun(1 + dv, 1 + du, dn, dm) * (double) (c == 0);
      }
    }
  double e1 = (double) EI1, e2 = (double) EI2;
  double V11 = e1 + (double) sAA - e1 * e1;
  double V22 = e2 + (double) sBB - e2 * e2;
  double V12 = (double) sAB - e1 * e2;
  NumericMatrix V(2, 2);
  V(0, 0) = V11; V(1, 0) = V12; V(0, 1) = V12; V(1, 1) = V22;
  return List::create(_["theta"] = NumericVector::create(e1, e2),
                      _["V"] = V);
}

// ── Exact attained p-values of S and T (Proposition S.2, dynamic programs) ──
typedef std::vector<long double> poly;

static poly pmul(const poly &a, const poly &b) {
  poly r(a.size() + b.size() - 1, 0.0L);
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] == 0) continue;
    for (size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
  }
  return r;
}
static poly padd(const poly &a, const poly &b) {
  poly r(std::max(a.size(), b.size()), 0.0L);
  for (size_t i = 0; i < a.size(); ++i) r[i] += a[i];
  for (size_t i = 0; i < b.size(); ++i) r[i] += b[i];
  return r;
}
static inline long double coef(const poly &p, int k) {
  return (k >= 0 && k < (int) p.size()) ? p[k] : 0.0L;
}

// [[Rcpp::export]]
NumericVector attained_pvalues(IntegerMatrix edges, int n, int m) {
  const int N = n + m;
  if (n < 1 || m < 1) stop("n and m must be positive");
  if (edges.nrow() != N - 1) stop("edges must form a spanning tree on n + m vertices");
  std::vector<std::vector<int> > nb(N);
  for (int e = 0; e < N - 1; ++e) {
    int a = edges(e, 0) - 1, c = edges(e, 1) - 1;
    if (a < 0 || a >= N || c < 0 || c >= N) stop("edge index out of range");
    nb[a].push_back(c); nb[c].push_back(a);
  }
  // root at vertex 0; order in which every vertex precedes its children
  std::vector<int> order, parent(N, -2), size(N, 1);
  order.reserve(N);
  parent[0] = -1; order.push_back(0);
  for (size_t i = 0; i < order.size(); ++i)
    for (int w : nb[order[i]])
      if (parent[w] == -2) { parent[w] = order[i]; order.push_back(w); }
  if ((int) order.size() != N) stop("edges do not form a spanning tree");

  const poly one{1.0L}, zpoly{0.0L, 1.0L}, zero{0.0L};
  std::vector<poly> f(N, zpoly), a0(N, zpoly), a1(N, zero), a2(N, one),
    D10(N, zpoly), D11(N, zero), D20(N, one), D21(N, zero);
  poly c(N + 1, 0.0L);
  long double s = 0;
  // children before parents: process order in reverse; when v is reached all
  // its children have been merged into it
  for (int i = N - 1; i >= 0; --i) {
    int v = order[i];
    for (int k = 0; k < (int) f[v].size(); ++k) c[k] += f[v][k];
    if (v == 0) break;
    int p = parent[v];
    if (size[v] == n || size[v] == m) s += 1;
    size[p] += size[v];
    // connected subsets: f_p <- f_p (1 + f_v)
    f[p] = pmul(f[p], padd(one, f[v]));
    // subsets without isolated vertices
    poly S1 = padd(padd(a0[v], a1[v]), a2[v]), A01 = padd(a0[v], a1[v]);
    poly na0 = pmul(a0[p], a2[v]);
    poly na1 = padd(pmul(a1[p], S1), pmul(a0[p], A01));
    poly na2 = pmul(a2[p], padd(a1[v], a2[v]));
    a0[p].swap(na0); a1[p].swap(na1); a2[p].swap(na2);
    // both samples without isolated vertices (label 1 marked by z)
    poly c1 = padd(D10[v], D11[v]), c2 = padd(D20[v], D21[v]);
    poly n10 = pmul(D10[p], D21[v]);
    poly n11 = padd(pmul(D11[p], D21[v]), pmul(padd(D10[p], D11[p]), c1));
    poly n20 = pmul(D20[p], D11[v]);
    poly n21 = padd(pmul(D21[p], D11[v]), pmul(padd(D20[p], D21[p]), c2));
    D10[p].swap(n10); D11[p].swap(n11); D20[p].swap(n20); D21[p].swap(n21);
    poly().swap(f[v]); poly().swap(a0[v]); poly().swap(a1[v]); poly().swap(a2[v]);
    poly().swap(D10[v]); poly().swap(D11[v]); poly().swap(D20[v]); poly().swap(D21[v]);
  }
  poly u = padd(a1[0], a2[0]), ub = padd(D11[0], D21[0]);
  long double binom = 1;
  for (int k = 1; k <= n; ++k) binom = binom * (N - n + k) / k;
  long double piS = (coef(c, n) + coef(c, m) - (n == m ? 2 : 1) * s) / binom;
  long double piT = (coef(u, n) + coef(u, m) - coef(ub, n)) / binom;
  return NumericVector::create(_["pi_S"] = (double) piS, _["pi_T"] = (double) piT);
}
