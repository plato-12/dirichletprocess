# Save this as profiling/gaussian_profiling.R
library(dirichletprocess)
library(profvis)

# Basic profiling
profile_basic <- function(n_obs = 500, n_iter = 100) {
  set.seed(123)
  y <- rnorm(n_obs)

  profvis({
    dp <- DirichletProcessGaussian(y)
    dp <- Fit(dp, n_iter)
  })
}

# More complex profiling
profile_complex <- function(n_obs = 300, n_iter = 100) {
  set.seed(456)
  y <- c(rnorm(n_obs/2, -2, 0.5), rnorm(n_obs/2, 2, 0.5))

  profvis({
    dp <- DirichletProcessGaussian(y)
    dp <- Fit(dp, n_iter)
  })
}

# Multivariate normal profiling
profile_mvn <- function(n_obs = 200, n_iter = 50) {
  set.seed(789)
  require(mvtnorm)
  y1 <- rmvnorm(n_obs/2, c(-2, -2), diag(2))
  y2 <- rmvnorm(n_obs/2, c(2, 2), diag(2))
  y <- rbind(y1, y2)

  profvis({
    dp <- DirichletProcessMvnormal(y)
    dp <- Fit(dp, n_iter)
  })
}

profile_basic()

# Run with: profile_basic(), profile_complex(), or profile_mvn()
