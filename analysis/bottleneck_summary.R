# R Script to Generate Profiling Results and Bottleneck Analysis Report as Markdown

analyze_bottlenecks <- function(output_filepath = "analysis/profiling_results/bottleneck_analysis.md") {

  # Ensure the output directory exists
  output_dir <- dirname(output_filepath)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created directory:", output_dir, "\n")
  }

  # Open a connection to the output file
  sink(output_filepath)

  # --- Start of Markdown Content ---

  cat("# dirichletprocess Package: Profiling Results and Bottleneck Analysis\n\n")
  cat("**Author**: Analysis based on profiling scripts\n")
  # To make the date dynamic if this script is run, otherwise use a static one
  cat("**Date**: ", format(Sys.Date(), "%B %d, %Y"), "\n\n", sep="") # Or use a static date like "May 19, 2025"

  cat("## 1. Introduction\n\n")
  cat("This report summarizes the performance profiling results for the `dirichletprocess` R package. The goal of this profiling was to identify computational bottlenecks within the core MCMC sampling algorithms, specifically to guide the C++ reimplementation efforts. Profiling was conducted using `profvis` for overall function tracing and `microbenchmark` for timing specific components. Different mixture types (Gaussian, Beta, MVN) and varying data sizes were considered.\n\n")
  cat("The key outputs analyzed are:\n\n")
  cat("* `benchmark_summary.txt`: Contains `microbenchmark` timings for `ClusterComponentUpdate`, `ClusterParameterUpdate`, and `UpdateAlpha`.\n\n")
  cat("* `profvis` HTML outputs (e.g., `gaussian_n100.html`, `beta_n100.html`, etc.): Provide detailed flame graphs and data views of the `Fit` function execution.\n\n")

  cat("## 2. Microbenchmark Component Analysis\n\n")
  cat("The `microbenchmark` results provide precise timings for the core iterative components of the MCMC sampler.\n\n")
  cat("**(Note: The following section reproduces the key insights from the provided `benchmark_summary.txt`)**\n\n")
  cat("```text\n") # Start of text code block
  cat("# Benchmark Results\n\n\n")
  cat("##  Gaussian ClusterComponentUpdate (n=200) \n")
  cat("                                       expr    min     lq     mean\n")
  cat("1 ClusterComponentUpdate(dp_gaussian_bench) 9.0073 9.1943 9.592545\n")
  cat("   median     uq     max neval\n")
  cat("1 9.33175 9.4064 14.6837    20\n\n\n")
  cat("##  Gaussian ClusterParameterUpdate (n=200) \n")
  cat("                                       expr   min     lq    mean\n")
  cat("1 ClusterParameterUpdate(dp_gaussian_bench) 106.1 107.75 117.405\n")
  cat("  median    uq   max neval\n")
  cat("1 108.65 113.1 255.4    20\n\n\n")
  cat("##  Gaussian UpdateAlpha (n=200) \n")
  cat("                            expr  min   lq  mean median   uq   max\n")
  cat("1 UpdateAlpha(dp_gaussian_bench) 16.5 16.7 19.76     17 17.4 106.7\n")
  cat("  neval\n")
  cat("1    50\n\n\n")
  cat("##  Beta ClusterComponentUpdate (n=100, Non-conjugate) \n")
  cat("                                   expr    min     lq    mean\n")
  cat("1 ClusterComponentUpdate(dp_beta_bench) 7.2424 7.2891 7.41192\n")
  cat("  median     uq    max neval\n")
  cat("1 7.3018 7.4016 8.1218    10\n\n\n")
  cat("##  Beta ClusterParameterUpdate (n=100, Non-conjugate) \n")
  cat("                                   expr    min     lq     mean\n")
  cat("1 ClusterParameterUpdate(dp_beta_bench) 488.15 489.75 498.9444\n")
  cat("  median     uq    max neval\n")
  cat("1 491.40 509.00 516.21     5\n\n\n")
  cat("##  Beta UpdateAlpha (n=100, Non-conjugate) \n")
  cat("                            expr    min     lq   mean median    uq\n")
  cat("1 UpdateAlpha(dp_beta_bench) 18.1 18.35 20.494  18.85 19.55\n")
  cat("    max neval\n")
  cat("1 113.4    50\n")
  cat("```\n\n") # End of text code block

  cat("### Key Findings from Microbenchmarks:\n\n")
  cat("#### Gaussian Mixture (Conjugate, n=200)\n\n")
  cat("* **`UpdateAlpha`**: Median ~17.0 µs. Very fast, as expected.\n\n")
  cat("* **`ClusterComponentUpdate`**: Median ~9.33 ms. This is a significant portion of the work, involving iteration over data points and clusters.\n\n")
  cat("* **`ClusterParameterUpdate`**: Median ~108.65 ms. This is the most time-consuming component for the conjugate Gaussian model, involving posterior draws for each active cluster.\n\n")
  cat("#### Beta Mixture (Non-conjugate, n=100)\n\n")
  cat("* **`UpdateAlpha`**: Median ~18.85 µs. Remains very fast.\n\n")
  cat("* **`ClusterComponentUpdate`**: Median ~7.30 ms. Comparable in impact to the Gaussian case, considering the smaller `n`.\n\n")
  cat("* **`ClusterParameterUpdate`**: Median ~491.40 ms. This is **extremely slow**, nearly 0.5 seconds per call for only 100 data points. It's roughly 4.5 times slower than its conjugate Gaussian counterpart, even with half the data size.\n\n")
  cat("**Conclusion from Microbenchmarks:**\n")
  cat("The `ClusterParameterUpdate` step is the most computationally intensive, especially for non-conjugate mixtures like Beta, where it's an order of magnitude slower. `ClusterComponentUpdate` is the second major contributor to runtime. `UpdateAlpha` is negligible in comparison.\n\n")

  cat("## 3. `profvis` Full Sampler Analysis (Inferred)\n\n")
  cat("The `profvis` HTML outputs (e.g., `gaussian_n1000.html`, `beta_n500.html`) provide a holistic view of the `Fit` function. While not directly embedded here, common patterns emerge:\n\n")
  cat("* The vast majority of the execution time within `Fit()` is spent inside its main MCMC loop.\n\n")
  cat("* Within this loop, the calls to `ClusterParameterUpdate()` and `ClusterComponentUpdate()` (and their S3 dispatched methods) will dominate the flame graph.\n\n")
  cat("### Conjugate Models (Gaussian, MVN)\n\n")
  cat("* **`ClusterParameterUpdate.*`**:\n\n")
  cat("  * `PosteriorDraw.*`: Calculating posterior hyperparameters (which involves iterating/summing over data points within each cluster) and then sampling from the known posterior distribution. For MVN models, matrix operations within `PosteriorDraw.mvnormal` (e.g., sums of squares and cross-products, inversions if not careful, Cholesky decompositions) would be significant.\n\n")
  cat("* **`ClusterComponentUpdate.*`**:\n\n")
  cat("  * `Likelihood.*`: These functions are called $N \\times K_{avg}$ times per MCMC iteration (where $N$ is the number of data points, $K_{avg}$ is the average number of clusters). For MVN, `mvtnorm::dmvnorm` is a key function.\n\n") # Escaped backslash for LaTeX
  cat("  * `Predictive.*`: Called $N$ times per MCMC iteration for evaluating new cluster probabilities.\n\n")
  cat("  * `ClusterLabelChange.*`: Management of R lists for `clusterParameters` and `pointsPerCluster` (especially resizing and subsetting) can introduce overhead, particularly if the number of clusters changes frequently.\n\n")
  cat("### Non-Conjugate Models (Beta)\n\n")
  cat("* **`ClusterParameterUpdate.nonconjugate`**: This is the **overwhelming bottleneck**.\n\n")
  cat("  * The primary reason is the call to `PosteriorDraw.nonconjugate`, which internally uses `MetropolisHastings`.\n\n")
  cat("  * **`MetropolisHastings` Loop**: This loop runs for `dpObj$mhDraws` (e.g., 250) iterations *for each active cluster* in *each MCMC iteration of `Fit`*.\n\n")
  cat("    * **`Likelihood.beta`**: Critically, this is called for *all data points assigned to the current cluster* at *every single Metropolis-Hastings step*. Summing log-likelihoods across these points is the most intensive part.\n\n")
  cat("    * `PriorDensity.beta`: Called once per MH step.\n\n")
  cat("    * `MhParameterProposal.beta`: Called once per MH step.\n\n")
  cat("  * **Low Metropolis-Hastings Acceptance Ratios**: The console output showed acceptance ratios around 0.3% - 0.9% for Beta mixtures. This means over 99% of the computationally expensive MH proposals (including likelihood calculations) are rejected. This highlights an inefficiency in the MH sampler's proposal/tuning for the Beta model, making the already expensive parameter update step even slower in terms of effective samples.\n\n")
  cat("* **`ClusterComponentUpdate.nonconjugate`**:\n\n")
  cat("  * `Likelihood.beta`: Still called frequently ($N \\times K_{avg}$ for existing clusters, and $N \\times m_{aux}$ for auxiliary parameters).\n\n") # Escaped backslash
  cat("  * `PriorDraw.beta`: Used for drawing auxiliary parameters.\n\n")
  cat("  * `ClusterLabelChange.nonconjugate`.\n\n")

  cat("## 4. Summary of Identified Bottlenecks for C++ Optimization\n\n")
  cat("Based on the microbenchmarks and the expected behavior from `profvis` analysis, the following are the primary targets for C++ optimization, in approximate order of impact:\n\n")
  cat("1. **`ClusterParameterUpdate` (Non-Conjugate Models)**:\n\n")
  cat("   * The entire **`MetropolisHastings` sampling loop** within `PosteriorDraw.nonconjugate`.\n\n")
  cat("   * Specifically, the repeated calculation of **`Likelihood.*`** (e.g., `Likelihood.beta`) over all data points in a cluster, for each MH step.\n\n")
  cat("   * Also, `PriorDensity.*` and `MhParameterProposal.*` within the MH loop.\n\n")
  cat("   * *Goal*: Make each MH step significantly faster.\n\n")
  cat("2. **`Likelihood.*` Functions (All Kernels)**:\n\n")
  cat("   * These are the most frequently called low-level computational routines.\n\n")
  cat("   * Optimizing these in C++ will benefit both `ClusterComponentUpdate` and `ClusterParameterUpdate` (especially non-conjugate).\n\n")
  cat("   * Focus on efficient mathematical computations and vectorized operations where possible.\n\n")
  cat("3. **`ClusterComponentUpdate` Core Logic (All Models)**:\n\n")
  cat("   * The outer loop over $N$ data points.\n\n")
  cat("   * The inner loop over $K$ clusters for calculating probabilities (heavy use of `Likelihood`).\n\n")
  cat("   * Calculation of new cluster probabilities (using `Predictive` for conjugate, or `Likelihood` with auxiliary parameters for non-conjugate).\n\n")
  cat("   * The categorical sampling step.\n\n")
  cat("4. **`PosteriorDraw.*` and `Predictive.*` (Conjugate Models)**:\n\n")
  cat("   * The calculation of posterior hyperparameters (often involves sums/means over data in clusters).\n\n")
  cat("   * The actual random sampling from these posteriors.\n\n")
  cat("   * Mathematical formulas for predictive distributions.\n\n")
  cat("   * For MVN, efficient matrix algebra is key.\n\n")
  cat("5. **`ClusterLabelChange.*`**:\n\n")
  cat("   * Efficient management of R data structures (cluster parameters, labels, counts) when clusters are added, removed, or re-indexed. C++ data structures like `std::vector` can offer better performance for dynamic resizing and access compared to R lists if large parts of this logic are moved to C++.\n\n")
  cat("6. **Data Structure Access and Management**:\n\n")
  cat("   * Minimizing copying of data between R and C++.\n\n")
  cat("   * Using efficient C++ data structures internally.\n\n")
  cat("   * Reducing overhead from repeated S3 method dispatch for core computational routines by having C++ handle larger chunks of work.\n\n")

  cat("## 5. Conclusion\n\n")
  cat("The profiling results clearly indicate that the `ClusterParameterUpdate` step, particularly for non-conjugate mixtures due to the nested Metropolis-Hastings sampling, is the most significant performance bottleneck. The `Likelihood` functions, being called extensively within both `ClusterParameterUpdate` and `ClusterComponentUpdate`, are also critical targets.\n\n")
  cat("Reimplementing these core components and their dependencies (like likelihood, prior, and posterior draw functions) in C++ is expected to yield substantial performance improvements. Additionally, addressing the low MH acceptance ratio for the Beta mixture through better proposal tuning (an R-level concern) would complement the C++ speedups for practical usability.\n")

  # --- End of Markdown Content ---
  sink() # Close the file connection
  cat("Markdown report generated at:", output_filepath, "\n")
}

# Example usage:
# Ensure the 'analysis/profiling_results/' directory exists or will be created by the function.
# analyze_bottlenecks()
# To run it, uncomment the line above and ensure your working directory is appropriate,
# or provide a full path for output_filepath.
# For example:
# analyze_bottlenecks(output_filepath = "my_profiling_report.md")

cat("R script 'analyze_bottlenecks' is defined. Call it to generate the Markdown file.\n")
analyze_bottlenecks()
