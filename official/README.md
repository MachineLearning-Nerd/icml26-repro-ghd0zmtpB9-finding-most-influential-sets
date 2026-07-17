# Finding Most Influential Sets

R implementation of algorithms for most influential set selection (MISS) in regression models using the Dinkelbach algorithm.

## Setup

```r
# Source functions
for (f in list.files("R", "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
# Optional Dependencies
install.packages("Rcpp") # Fast Top-K
install.packages(c("gbm", "xgboost", "randomForest")) # PLM options (gradient boosting, random forest)
devtools::install_github("nk027/infuential_sets") # Greedy algorithm
```

## Usage

```r
# Obtain data
N <- 1000L
x <- rnorm(N)
y <- x + rnorm(N)
# Create model
model <- lm(y ~ 0 + x)

# Find most influential set of size k
result <- find_miss(model, k = 10)
# Find sets for sizes 1 to K
results <- find_misses(model, K = 5)
```

### Available Functions

- `find_miss(model, k)`: Find optimal set of size k
- `enumerate_miss(model, k)`: Exhaustive enumeration
- `greedy_miss(model, k)`: Greedy approximation

