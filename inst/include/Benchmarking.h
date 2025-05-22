#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include <RcppArmadillo.h>
#include <chrono>
#include <vector>
#include <string>
#include "DirichletProcess.h"

namespace dp {

// Timer class for C++ benchmarking
class Timer {
private:
  std::chrono::high_resolution_clock::time_point start_time;

public:
  // Start timer
  void start() {
    start_time = std::chrono::high_resolution_clock::now();
  }

  // Return elapsed time in milliseconds
  double elapsed_ms() {
    auto end_time = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(
        end_time - start_time).count();
  }
};

// Memory tracker class definition
class MemoryTracker {
private:
  std::vector<size_t> allocations;
  std::vector<std::string> descriptions;

public:
  // Record an allocation
  void record(size_t bytes, const std::string& description) {
    allocations.push_back(bytes);
    descriptions.push_back(description);
  }

  // Get total allocated memory
  size_t total() const {
    size_t total = 0;
    for (size_t alloc : allocations) {
      total += alloc;
    }
    return total;
  }

  // Get summary as R data frame
  Rcpp::DataFrame summary() const {
    if (allocations.empty()) {
      return Rcpp::DataFrame::create();
    }

    Rcpp::NumericVector bytes(allocations.begin(), allocations.end());
    Rcpp::CharacterVector desc(descriptions.begin(), descriptions.end());

    return Rcpp::DataFrame::create(
      Rcpp::Named("description") = desc,
      Rcpp::Named("bytes") = bytes,
      Rcpp::Named("mb") = bytes / (1024.0 * 1024.0)
    );
  }

  // Clear all records
  void clear() {
    allocations.clear();
    descriptions.clear();
  }
};

// Function to benchmark a component
template<typename F>
double benchmark_function(F func, int times = 10) {
  Timer timer;
  std::vector<double> timings(times);

  for (int i = 0; i < times; i++) {
    timer.start();
    func();
    timings[i] = timer.elapsed_ms();
  }

  // Calculate median
  std::sort(timings.begin(), timings.end());
  if (times % 2 == 0) {
    return (timings[times/2 - 1] + timings[times/2]) / 2.0;
  } else {
    return timings[times/2];
  }
}

// FUNCTION DECLARATIONS ONLY (not definitions)
// Forward declare functions that will be defined in the .cpp file
size_t current_memory_usage();
Rcpp::List benchmark_cpp_components(const Rcpp::List& dpObj,
                                    const Rcpp::StringVector& components,
                                    int times = 10);
Rcpp::DataFrame get_memory_tracking();
void clear_memory_tracking();

// Declare the global memory tracker (but don't define it here)
extern MemoryTracker g_memory_tracker;

} // namespace dp

#endif
