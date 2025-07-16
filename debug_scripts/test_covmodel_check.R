# Test covModel check
library(dirichletprocess)
set_use_cpp(FALSE)  # Turn off C++

cat("=== TESTING COVMODEL CHECK ===\n")

# Create E model
params_1d <- list(
  mu0 = 0,
  kappa0 = 1,
  nu = 3,
  Lambda = matrix(1, 1, 1),
  covModel = "E"
)

md <- MvnormalCreate(params_1d)
data_1d <- matrix(rnorm(10), ncol = 1)
dp <- DirichletProcessCreate(data_1d, md)

# Test the condition
cat("Testing covModel condition...\n")
cat(sprintf("exists('priorParameters', dp$mixingDistribution): %s\n", 
            exists("priorParameters", dp$mixingDistribution)))
cat(sprintf("!is.null(dp$mixingDistribution$priorParameters$covModel): %s\n", 
            !is.null(dp$mixingDistribution$priorParameters$covModel)))
cat(sprintf("dp$mixingDistribution$priorParameters$covModel: %s\n", 
            dp$mixingDistribution$priorParameters$covModel))

skip_expansion <- (exists("priorParameters", dp$mixingDistribution) && 
                   !is.null(dp$mixingDistribution$priorParameters$covModel) &&
                   dp$mixingDistribution$priorParameters$covModel %in% c("E", "V"))

cat(sprintf("skip_expansion: %s\n", skip_expansion))

# Full condition
will_expand <- (inherits(dp, "mvnormal") && !skip_expansion)
cat(sprintf("will_expand: %s\n", will_expand))

cat("\n=== COVMODEL CHECK TEST COMPLETE ===\n")