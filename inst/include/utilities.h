#ifndef DIRICHLETPROCESS_UTILITIES_H
#define DIRICHLETPROCESS_UTILITIES_H

#include <RcppArmadillo.h>

namespace dirichletprocess {

// Function to sample from a categorical distribution defined by a vector of probabilities.
// It returns the 0-based index of the chosen category.
inline int sample_categorical(const arma::vec& probs) {
  // Ensure probabilities sum to 1 (within a small tolerance)
  if (std::abs(arma::sum(probs) - 1.0) > 1e-8) {
    Rcpp::warning("Probabilities in sample_categorical do not sum to 1.");
  }

  // Generate a single uniform random number between 0 and 1
  double u = R::runif(0, 1);

  double cumulative_prob = 0.0;
  for (arma::uword i = 0; i < probs.n_elem; ++i) {
    cumulative_prob += probs[i];
    if (u < cumulative_prob) {
      return i; // Return the index of the first interval `u` falls into
    }
  }

  // Fallback for floating-point precision issues: return the last index.
  return probs.n_elem - 1;
}

} // namespace dirichletprocess

#endif // DIRICHLETPROCESS_UTILITIES_H
