# R/manual_mcmc_cpp.R

#' Manual C++ MCMC Runner Class
#'
#' @field ptr External pointer to C++ MCMCRunnerManual object
#' @field dp_obj Original Dirichlet process object
#' @field distribution_type Type of mixing distribution
#' @export
CppMCMCRunner <- setRefClass("CppMCMCRunner",
                             fields = list(
                               ptr = "externalptr",
                               dp_obj = "ANY",
                               distribution_type = "character"
                             ),

                             methods = list(
                               initialize = function(dp_object, n_iter = 1000, n_burn = 100, thin = 1) {
                                 # Prepare data
                                 data <- as.matrix(dp_object$data)

                                 # Prepare parameters (works for ALL distributions)
                                 mixing_params <- prepare_mixing_dist_params(dp_object)
                                 mcmc_params <- prepare_mcmc_params(dp_object, n_iter, TRUE, n_burn, thin)

                                 # Create C++ runner
                                 ptr <<- create_mcmc_runner_cpp(data, mixing_params, mcmc_params)
                                 dp_obj <<- dp_object
                                 distribution_type <<- mixing_params$type
                               },

                               step_assignments = function() {
                                 "Update cluster assignments using Algorithm 8"
                                 step_assignments_cpp(ptr)
                                 invisible(.self)
                               },

                               step_parameters = function() {
                                 "Update cluster parameters"
                                 step_parameters_cpp(ptr)
                                 invisible(.self)
                               },

                               step_concentration = function() {
                                 "Update concentration parameter"
                                 step_concentration_cpp(ptr)
                                 invisible(.self)
                               },

                               step = function() {
                                 "Perform one complete MCMC iteration"
                                 perform_iteration_cpp(ptr)
                                 invisible(.self)
                               },

                               get_state = function() {
                                 "Get current MCMC state"
                                 state <- get_state_cpp(ptr)

                                 # Convert parameters back to original format based on distribution
                                 if (distribution_type == "mvnormal") {
                                   # Special handling for multivariate normal
                                   state$cluster_params <- lapply(state$cluster_params, function(p) {
                                     list(mu = p[1:dp_obj$mixingDistribution$d],
                                          Sigma = matrix(p[-(1:dp_obj$mixingDistribution$d)],
                                                         nrow = dp_obj$mixingDistribution$d))
                                   })
                                 }

                                 state
                               },

                               get_results = function() {
                                 "Get accumulated MCMC results"
                                 results <- get_results_cpp(ptr)

                                 # Format results to match standard Fit() output
                                 results$numberClusters <- results$n_clusters
                                 results$clusterParameters <- results$cluster_params
                                 results$clusterLabels <- results$cluster_labels + 1  # R uses 1-based indexing
                                 results$alpha <- results$alpha
                                 results$data <- dp_obj$data
                                 results$mixingDistribution <- dp_obj$mixingDistribution

                                 class(results) <- c("dirichletprocess", "list")
                                 results
                               },

                               is_complete = function() {
                                 "Check if all iterations are complete"
                                 is_complete_cpp(ptr)
                               },

                               run_manual = function(callback = NULL, progress = TRUE,
                                                     custom_assignment_sampler = NULL,
                                                     custom_parameter_sampler = NULL) {
                                 "Run MCMC with optional callbacks and custom samplers"

                                 if (progress) {
                                   pb <- txtProgressBar(min = 0, max = 1, style = 3)
                                 }

                                 while (!is_complete()) {
                                   # Custom or standard cluster assignment
                                   if (!is.null(custom_assignment_sampler)) {
                                     state <- get_state()
                                     new_labels <- custom_assignment_sampler(state, dp_obj)
                                     set_labels_cpp(ptr, new_labels)
                                   } else {
                                     step_assignments()
                                   }

                                   # Custom or standard parameter update
                                   if (!is.null(custom_parameter_sampler)) {
                                     state <- get_state()
                                     new_params <- custom_parameter_sampler(state, dp_obj)
                                     # Would need set_params_cpp implementation
                                   } else {
                                     step_parameters()
                                   }

                                   # Always update concentration
                                   step_concentration()

                                   # User callback
                                   if (!is.null(callback)) {
                                     state <- get_state()
                                     callback(state, .self)
                                   }

                                   # Progress
                                   if (progress) {
                                     state <- get_state()
                                     setTxtProgressBar(pb, state$iteration / mcmc_params$n_iter)
                                   }
                                 }

                                 if (progress) close(pb)

                                 get_results()
                               }
                             )
)

#' Create Manual C++ MCMC Runner
#'
#' @param dp_obj Dirichlet process object (any distribution)
#' @param n_iter Number of iterations
#' @param n_burn Burn-in iterations
#' @param thin Thinning interval
#' @export
create_cpp_mcmc_runner <- function(dp_obj, n_iter = 1000, n_burn = 100, thin = 1) {

  # Check if C++ is available for this distribution
  if (!can_use_cpp(dp_obj)) {
    stop("C++ backend not available for distribution: ",
         class(dp_obj$mixingDistribution))
  }

  CppMCMCRunner$new(dp_obj, n_iter, n_burn, thin)
}
