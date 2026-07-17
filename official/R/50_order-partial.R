# Partial ordering: get top-k indices efficiently

# Use Rcpp if available, otherwise falls back to base R
if (requireNamespace("Rcpp", quietly = TRUE)) {
  Rcpp::sourceCpp(
    code = '
#include <Rcpp.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits>

using namespace Rcpp;

struct Node {
  double v;
  int idx;
};

// "Better" for decreasing: larger value first; tie: smaller idx first
struct BetterDec {
  inline bool operator()(const Node& a, const Node& b) const noexcept {
    if (a.v > b.v) return true;
    if (a.v < b.v) return false;
    return a.idx < b.idx;
  }
};

// "Better" for increasing: smaller value first; tie: smaller idx first
struct BetterInc {
  inline bool operator()(const Node& a, const Node& b) const noexcept {
    if (a.v < b.v) return true;
    if (a.v > b.v) return false;
    return a.idx < b.idx;
  }
};

template <class Better>
static inline IntegerVector topk_heap(const double* xp, R_xlen_t n, int k, const Better& better) {
  if (n == 0 || k <= 0) return IntegerVector();
  if (k > (int)n) k = (int)n;

  // k==1 fast path
  if (k == 1) {
    R_xlen_t best = 0;
    for (R_xlen_t i = 1; i < n; ++i) {
      Node a{xp[i], (int)i};
      Node b{xp[best], (int)best};
      if (better(a, b)) best = i;
    }
    IntegerVector out(1);
    out[0] = (int)best + 1;
    return out;
  }

  // priority_queue top is the element that is NOT "better" than any other => the WORST
  std::priority_queue<Node, std::vector<Node>, Better> pq(better);

  // fill first k
  for (int i = 0; i < k; ++i) pq.push(Node{xp[i], i});

  // scan remaining
  for (R_xlen_t i = (R_xlen_t)k; i < n; ++i) {
    Node cand{xp[i], (int)i};
    const Node& worst = pq.top();
    if (better(cand, worst)) {
      pq.pop();
      pq.push(cand);
    }
  }

  // extract and sort best->worst
  std::vector<Node> res;
  res.reserve((size_t)k);
  while (!pq.empty()) {
    res.push_back(pq.top());
    pq.pop();
  }

  std::sort(res.begin(), res.end(),
            [&](const Node& a, const Node& b) { return better(a, b); });

  IntegerVector out(k);
  for (int i = 0; i < k; ++i) out[i] = res[(size_t)i].idx + 1;
  return out;
}

// [[Rcpp::export]]
IntegerVector order_partial_cpp(SEXP xSEXP, int k, bool decreasing) {
  if (TYPEOF(xSEXP) != REALSXP) stop("Input must be a double (REALSXP)");
  const R_xlen_t n = XLENGTH(xSEXP);
  if (n > (R_xlen_t)std::numeric_limits<int>::max())
    stop("x is too long: indices must fit in 32-bit integer");

  const double* xp = REAL(xSEXP);

  if (decreasing) {
    return topk_heap(xp, n, k, BetterDec{});
  } else {
    return topk_heap(xp, n, k, BetterInc{});
  }
}
'
  )
  order_partial <- function(x, k, decreasing = FALSE) {
    order_partial_cpp(x, as.integer(k), isTRUE(decreasing))
  }
} else {
  message("Rcpp not available; using plain R for (partial) order.")
  order_partial <- function(x, k, decreasing = FALSE) {
    order(x, decreasing = decreasing)[seq_len(min(k, length(x)))]
  }
}
