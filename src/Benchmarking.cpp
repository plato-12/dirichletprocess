// src/Benchmarking.cpp
#include "../inst/include/Benchmarking.h"
#include "../inst/include/RcppConversions.h"
#include "../inst/include/NormalDistribution.h"
#include "../inst/include/BetaDistribution.h"
#include "../inst/include/MVNormalDistribution.h"
#include "../inst/include/MVNormal2Distribution.h"
#include "../inst/include/WeibullDistribution.h"
#include "../inst/include/ExponentialDistribution.h"
#include "../inst/include/HierarchicalDP.h"
#include "../inst/include/MarkovDP.h"

// Platform-specific includes for memory tracking
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

#ifdef __APPLE__
#include <mach/mach.h>
#include <sys/resource.h>
#endif

#ifdef __linux__
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace dp {

// Global memory tracker instance
MemoryTracker g_memory_tracker;

// [[Rcpp::export]]
size_t current_memory_usage() {
#ifdef _WIN32
  // Windows implementation
  PROCESS_MEMORY_COUNTERS_EX pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
    return pmc.WorkingSetSize;
  }
  return 0;
#elif defined(__APPLE__) || defined(__linux__)
  // Mac/Linux implementation
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) == 0) {
#ifdef __APPLE__
    return usage.ru_maxrss; // macOS gives bytes
#else
    return usage.ru_maxrss * 1024; // Linux gives KB
#endif
  }
  return 0;
#else
  // Fallback for unsupported platforms
  return 0;
#endif
}

// [[Rcpp::export]]
Rcpp::List benchmark_cpp_components(const Rcpp::List& dpObj,
                                    const Rcpp::StringVector& components,
                                    int times) {
  Rcpp::List results;

  try {
    // Convert R object to C++
    DirichletProcess* dp = createDPFromR(dpObj);

    // Run benchmarks for each requested component
    for (int i = 0; i < components.size(); i++) {
      std::string component = Rcpp::as<std::string>(components[i]);
      double time_ms = 0.0;

      if (component == "clusterComponentUpdate") {
        time_ms = benchmark_function([&]() { dp->clusterComponentUpdate(); }, times);
      } else if (component == "clusterParameterUpdate") {
        time_ms = benchmark_function([&]() { dp->clusterParameterUpdate(); }, times);
      } else if (component == "updateAlpha") {
        time_ms = benchmark_function([&]() { dp->updateAlpha(); }, times);
      } else if (component == "likelihood") {
        // We need to handle likelihood benchmarking separately as it depends on the distribution type
        std::string distType = getDistributionType(dpObj);
        Rcpp::List mdObj = dpObj["mixingDistribution"];

        if (distType == "normal") {
          // Fix: Properly convert Rcpp types
          Rcpp::NumericVector priorParams = mdObj["priorParameters"];
          NormalMixingDistribution* md = new NormalMixingDistribution(priorParams);
          arma::vec x = convertVector(Rcpp::as<Rcpp::NumericVector>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() { md->likelihood(x, theta); }, times);

          delete md;
        } else if (distType == "beta") {
          Rcpp::NumericVector priorParams = mdObj["priorParameters"];
          BetaMixingDistribution* md = new BetaMixingDistribution(priorParams);
          arma::vec x = convertVector(Rcpp::as<Rcpp::NumericVector>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() { md->likelihood(x, theta); }, times);

          delete md;
        } else if (distType == "mvnormal") {
          Rcpp::List priorParams = mdObj["priorParameters"];
          MVNormalMixingDistribution* md = new MVNormalMixingDistribution(priorParams);
          arma::mat x = convertMatrix(Rcpp::as<Rcpp::NumericMatrix>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() {
            md->likelihood(arma::vectorise(x), theta);
          }, times);

          delete md;
        } else if (distType == "mvnormal2") {
          Rcpp::List priorParams = mdObj["priorParameters"];
          MVNormal2MixingDistribution* md = new MVNormal2MixingDistribution(priorParams);
          arma::mat x = convertMatrix(Rcpp::as<Rcpp::NumericMatrix>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() {
            md->likelihood(arma::vectorise(x), theta);
          }, times);

          delete md;
        } else if (distType == "weibull") {
          Rcpp::NumericVector priorParams = mdObj["priorParameters"];
          Rcpp::NumericVector mhStep = mdObj["mhStepSize"];
          Rcpp::NumericVector hyperPrior = mdObj.containsElementNamed("hyperPriorParameters") ?
          Rcpp::as<Rcpp::NumericVector>(mdObj["hyperPriorParameters"]) :
            Rcpp::NumericVector::create();

          WeibullMixingDistribution* md = new WeibullMixingDistribution(priorParams, mhStep, hyperPrior);
          arma::vec x = convertVector(Rcpp::as<Rcpp::NumericVector>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() { md->likelihood(x, theta); }, times);

          delete md;
        } else if (distType == "exponential") {
          Rcpp::NumericVector priorParams = mdObj["priorParameters"];
          ExponentialMixingDistribution* md = new ExponentialMixingDistribution(priorParams);
          arma::vec x = convertVector(Rcpp::as<Rcpp::NumericVector>(dpObj["data"]));
          Rcpp::List theta = dpObj["clusterParameters"];

          time_ms = benchmark_function([&]() { md->likelihood(x, theta); }, times);

          delete md;
        } else {
          Rcpp::warning("Unsupported distribution type for likelihood benchmarking: " + distType);
        }
      } else {
        Rcpp::warning("Unknown component: " + component);
        continue;
      }

      results[component] = time_ms;
    }

    // Clean up
    delete dp;

  } catch (std::exception& e) {
    Rcpp::stop("Error in C++ benchmarking: " + std::string(e.what()));
  }

  return results;
}

// [[Rcpp::export]]
Rcpp::DataFrame get_memory_tracking() {
  return g_memory_tracker.summary();
}

// [[Rcpp::export]]
void clear_memory_tracking() {
  g_memory_tracker.clear();
}

} // namespace dp
