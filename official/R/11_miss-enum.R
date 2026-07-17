# Enumerate all possible sets
enumerate_miss <- function(model, k, sign = 1L, verbose = TRUE) {
  xr <- get_lm_xr(model)
  X <- xr$X * sign
  R <- xr$R

  if (isTRUE(verbose)) {
    message(
      "Runtime (at 1 us per iteration): ",
      to_time(choose(NROW(X), k), unit = "us")
    )
  }

  num <- X * R
  den <- X^2
  combs <- combn(length(num), k, simplify = FALSE)
  vals <- vapply(
    combs,
    \(set, A, B, C) {
      sum(A[set]) / (C - sum(B[set]))
    },
    numeric(1L),
    A = num,
    B = den,
    C = sum(den)
  )
  best_idx <- which.max(vals)
  list(
    k = k,
    best_S = combs[[best_idx]],
    best_value = coef(model) - vals[best_idx] * sign
  )
}
