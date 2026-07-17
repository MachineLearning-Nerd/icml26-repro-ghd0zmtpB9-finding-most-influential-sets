#!/usr/bin/env Rscript

# Executes the released implementation unchanged after sourcing its R files.
root <- normalizePath(".")
out_dir <- file.path(root, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in list.files(file.path(root, "official", "R"), "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

set.seed(260605919L)

as_set <- function(x) paste(sort(as.integer(x)), collapse = ";")

case_rows <- list()
data_rows <- list()
case_id <- 0L
for (n in c(8L, 10L, 12L, 14L, 16L, 18L)) {
  for (k in seq_len(min(5L, n - 2L))) {
    for (rep in seq_len(4L)) {
      case_id <- case_id + 1L
      x <- rnorm(n)
      y <- 0.3 * x + (1 + 0.2 * abs(x)) * rt(n, df = 4L)
      model <- lm(y ~ 0 + x)
      official <- find_miss(model, k = k)
      enumerated <- enumerate_miss(model, k = k, verbose = FALSE)
      case_rows[[case_id]] <- data.frame(
        case_id = case_id, n = n, k = k, rep = rep,
        official_set = as_set(official$best_S), official_value = as.numeric(official$best_value),
        enumerated_set = as_set(enumerated$best_S), enumerated_value = as.numeric(enumerated$best_value),
        iterations = official$iter
      )
      data_rows[[case_id]] <- data.frame(case_id = case_id, index = seq_len(n), x = x, y = y)
    }
  }
}
write.csv(do.call(rbind, case_rows), file.path(out_dir, "official_cases.csv"), row.names = FALSE)
write.csv(do.call(rbind, data_rows), file.path(out_dir, "official_case_data.csv"), row.names = FALSE)

scale_rows <- list()
for (i in seq_along(c(1000L, 10000L, 100000L, 1000000L))) {
  n <- c(1000L, 10000L, 100000L, 1000000L)[i]
  x <- rnorm(n)
  y <- x + rnorm(n)
  model <- lm(y ~ 0 + x)
  elapsed <- system.time(for (repeat_id in seq_len(5L)) result <- find_miss(model, k = min(100L, n - 2L)))[["elapsed"]] / 5L
  scale_rows[[i]] <- data.frame(n = n, k = min(100L, n - 2L), seconds = elapsed, repetitions = 5L, iterations = result$iter,
                                 best_value = as.numeric(result$best_value))
}
write.csv(do.call(rbind, scale_rows), file.path(out_dir, "official_scaling.csv"), row.names = FALSE)

wine <- read.csv(file.path(root, "official", "data", "winequality-red.csv"))
wine_model <- lm(quality ~ 0 + alcohol, data = wine)
wine_result <- find_miss(wine_model, k = 25L)
write.csv(data.frame(n = nrow(wine), k = 25L, iterations = wine_result$iter,
                     best_value = as.numeric(wine_result$best_value), best_set = as_set(wine_result$best_S)),
          file.path(out_dir, "official_realdata.csv"), row.names = FALSE)

# Theorem 2 setup: oracle residualized inputs plus a controlled estimated-nuisance perturbation.
sep_rows <- list()
row_id <- 0L
for (rep in seq_len(25L)) {
  n <- 20L
  x_oracle <- rnorm(n)
  # Strong but finite outliers make a unique, measurable oracle separation gap likely.
  y_oracle <- 0.6 * x_oracle + rt(n, df = 5L)
  y_oracle[c(n - 3L, n - 1L)] <- y_oracle[c(n - 3L, n - 1L)] + c(6, -5)
  oracle_model <- lm(y_oracle ~ 0 + x_oracle)
  oracle <- find_miss(oracle_model, k = 4L)
  xr <- get_lm_xr(oracle_model)
  w <- xr$X * xr$R
  cst <- xr$X^2
  total <- sum(cst)
  all_sets <- combn(n, 4L, simplify = FALSE)
  all_values <- vapply(all_sets, function(S) sum(w[S]) / (total - sum(cst[S])), numeric(1L))
  ordered <- sort(all_values, decreasing = TRUE)
  gap <- ordered[1L] - ordered[2L]
  dx <- as.vector(scale(sin(seq_len(n) * (rep + 1L))))
  dy <- as.vector(scale(cos(seq_len(n) * (rep + 3L))))
  for (alpha in c(0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2)) {
    row_id <- row_id + 1L
    estimated_model <- lm((y_oracle + alpha * dy) ~ 0 + I(x_oracle + alpha * dx))
    estimated <- find_miss(estimated_model, k = 4L)
    sep_rows[[row_id]] <- data.frame(
      rep = rep, alpha = alpha, oracle_set = as_set(oracle$best_S), estimated_set = as_set(estimated$best_S),
      exact_recovery = identical(sort(oracle$best_S), sort(estimated$best_S)),
      oracle_value = as.numeric(oracle$best_value), estimated_value = as.numeric(estimated$best_value),
      separation_gap = gap, iterations = estimated$iter
    )
  }
}
write.csv(do.call(rbind, sep_rows), file.path(out_dir, "official_separation.csv"), row.names = FALSE)

# Negative control: rank-deficient deletion is rejected by the released positive-denominator guard.
control <- tryCatch({
  x <- c(1, 0, 0, 0)
  y <- c(1, 0, 0, 0)
  find_miss(lm(y ~ 0 + x), k = 1L)
  "UNEXPECTED_ACCEPT"
}, error = function(e) "REJECTED_NONPOSITIVE_DENOMINATOR")
writeLines(control, file.path(out_dir, "official_control.txt"))

cat("Official findingMIS reproduction complete\n")
