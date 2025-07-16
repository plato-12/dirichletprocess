# Debug E model dimensions
library(dirichletprocess)
set_use_cpp(TRUE)

cat("=== DEBUGGING E MODEL DIMENSIONS ===\n")

# Test parameter creation for E model
params_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1, 1, 1),
  covModel = "E"
)

# Create mixing distribution
md <- MvnormalCreate(params_1d)

# Test PriorDraw
cat("Testing PriorDraw...\n")
prior_draw <- PriorDraw(md, 2)
cat("PriorDraw result:\n")
print(prior_draw)

cat("\nDimensions of mu:\n")
print(dim(prior_draw$mu))
cat("Dimensions of sig:\n")
print(dim(prior_draw$sig))

# Test PosteriorDraw
cat("\nTesting PosteriorDraw...\n")
data_1d <- matrix(rnorm(10), ncol = 1)
posterior_draw <- PosteriorDraw(md, data_1d, 1)
cat("PosteriorDraw result:\n")
print(posterior_draw)

cat("\nDimensions of mu:\n")
print(dim(posterior_draw$mu))
cat("Dimensions of sig:\n")
print(dim(posterior_draw$sig))

cat("\n=== E MODEL DEBUG COMPLETE ===\n")