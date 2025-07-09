#' Mixing Distribution Likelihood
#'
#' Evaluate the Likelihood of some data \eqn{x} for some parameter \eqn{\theta}.
#'
#' @param mdObj Mixing Distribution
#' @param x Data
#' @param theta Parameters of distribution
#' @return Likelihood of the data
#' @export
Likelihood <- function(mdObj, x, theta) {
  if (using_cpp()) {
    # Get the distribution type
    dist_type <- class(mdObj)[class(mdObj) != "list" & class(mdObj) != "MixingDistribution"][1]

    tryCatch({
      # Convert x to numeric vector if needed
      x <- as.numeric(x)

      # Dispatch to the appropriate C++ function
      if (dist_type == "weibull") {
        # Extract parameter arrays
        alpha_array <- theta[[1]]
        lambda_array <- theta[[2]]

        # Get the number of clusters from the third dimension
        num_clusters <- dim(alpha_array)[3]
        if (is.null(num_clusters)) num_clusters <- 1

        # If there's only one cluster, call the C++ function once
        if (num_clusters == 1) {
          alpha <- as.numeric(alpha_array[1, 1, 1])
          lambda <- as.numeric(lambda_array[1, 1, 1])
          return(weibull_likelihood_cpp(x, alpha, lambda))
        } else {
          # Multiple clusters: return likelihood for each cluster
          result <- numeric(num_clusters)
          for (k in 1:num_clusters) {
            alpha <- as.numeric(alpha_array[1, 1, k])
            lambda <- as.numeric(lambda_array[1, 1, k])
            # For each cluster, calculate likelihood for all data points,
            # then take the first value (for single data point case)
            lik_values <- weibull_likelihood_cpp(x, alpha, lambda)
            result[k] <- lik_values[1]
          }
          return(result)
        }
      } else if (dist_type == "normal" || dist_type == "gaussian") {
        # Similar handling for normal distribution
        mu_array <- theta[[1]]
        sigma_array <- theta[[2]]

        num_clusters <- dim(mu_array)[3]
        if (is.null(num_clusters)) num_clusters <- 1

        if (num_clusters == 1) {
          mu <- as.numeric(mu_array[1, 1, 1])
          sigma <- as.numeric(sigma_array[1, 1, 1])
          return(normal_likelihood_cpp(x, mu, sigma))
        } else {
          result <- numeric(num_clusters)
          for (k in 1:num_clusters) {
            mu <- as.numeric(mu_array[1, 1, k])
            sigma <- as.numeric(sigma_array[1, 1, k])
            lik_values <- normal_likelihood_cpp(x, mu, sigma)
            result[k] <- lik_values[1]
          }
          return(result)
        }
      } else if (dist_type == "exponential") {
        # Handle exponential distribution
        lambda_array <- theta[[1]]

        num_clusters <- dim(lambda_array)[3]
        if (is.null(num_clusters)) num_clusters <- 1

        if (num_clusters == 1) {
          lambda <- as.numeric(lambda_array[1, 1, 1])
          return(exponential_likelihood_cpp(x, lambda))
        } else {
          result <- numeric(num_clusters)
          for (k in 1:num_clusters) {
            lambda <- as.numeric(lambda_array[1, 1, k])
            lik_values <- exponential_likelihood_cpp(x, lambda)
            result[k] <- lik_values[1]
          }
          return(result)
        }
      } else if (dist_type == "beta") {
        # Handle beta distribution
        alpha_array <- theta[[1]]
        beta_array <- theta[[2]]
        maxT <- ifelse(is.null(mdObj$maxT), 1, mdObj$maxT)

        num_clusters <- dim(alpha_array)[3]
        if (is.null(num_clusters)) num_clusters <- 1

        if (num_clusters == 1) {
          alpha <- as.numeric(alpha_array[1, 1, 1])
          beta <- as.numeric(beta_array[1, 1, 1])
          return(beta_likelihood_cpp(x, alpha, beta, maxT))
        } else {
          result <- numeric(num_clusters)
          for (k in 1:num_clusters) {
            alpha <- as.numeric(alpha_array[1, 1, k])
            beta <- as.numeric(beta_array[1, 1, k])
            lik_values <- beta_likelihood_cpp(x, alpha, beta, maxT)
            result[k] <- lik_values[1]
          }
          return(result)
        }
      } else {
        # For other distributions, fall back to R
        stop("C++ implementation not available for distribution: ", dist_type)
      }
    }, error = function(e) {
      warning("C++ implementation failed, falling back to R: ", e$message)
    })
  }

  # Original implementation (falls back to this if C++ is not enabled or fails)
  UseMethod("Likelihood", mdObj)
}
