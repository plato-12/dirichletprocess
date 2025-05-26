#' Fit a Hidden Markov Dirichlet Process Model
#'
#' @param dpObj Initialised Dirichlet Process object
#' @param its Number of iterations to use
#' @param updatePrior Logical flag, defaults to \code{FALSE}. Set whether the parameters of the base measure are updated.
#' @param progressBar Logical flag indicating whether to display a progress bar.
#' @return A Dirichlet Process object with the fitted cluster parameters and states.
#' @export
Fit.markov <- function(dpObj, its, updatePrior=F, progressBar = F){

  dpObj <- fit_hmm(dpObj, its, progressBar)

  return(dpObj)
}

# Helper function to ensure proper initialization of HMM parameters
initialize_hmm_params <- function(dpObj) {
  mdobj <- dpObj$mixingDistribution

  # Ensure all theta_k entries have proper structure
  if (inherits(mdobj, "normal")) {
    for (k in seq_along(mdobj$theta_k)) {
      if (!is.list(mdobj$theta_k[[k]]) ||
          !all(c("mean", "sd") %in% names(mdobj$theta_k[[k]]))) {
        # Convert to proper structure
        if (is.numeric(mdobj$theta_k[[k]])) {
          mdobj$theta_k[[k]] <- list(
            mean = mdobj$theta_k[[k]],
            sd = ifelse(is.null(mdobj$tau), 1, sqrt(1/mdobj$tau))
          )
        }
      }
    }
    dpObj$mixingDistribution <- mdobj
  }

  return(dpObj)
}

fit_hmm <- function(dpObj, its, progressBar=F){

  if (progressBar){
    pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
  }

  alphaChain <- numeric(its)
  betaChain <- numeric(its)
  statesChain <- vector("list", its)
  paramChain <- vector("list", its)

  # Ensure initial parameters are correctly structured
  dpObj <- initialize_hmm_params(dpObj)

  for(i in seq_len(its)){

    alphaChain[i] <- dpObj$alpha
    betaChain[i] <- dpObj$beta
    statesChain[[i]] <- dpObj$states


    paramChain[[i]] <- dpObj$uniqueParams

    dpObj <- UpdateStates(dpObj)
    dpObj <- UpdateAlphaBeta(dpObj)
    dpObj <- param_update(dpObj)


    if (progressBar) {
      setTxtProgressBar(pb, i)
    }

  }

  dpObj$alphaChain <- alphaChain
  dpObj$betaChain <- betaChain
  dpObj$statesChain <- statesChain
  dpObj$paramChain <- paramChain
  if (progressBar) {
    close(pb)
  }

  return(dpObj)
}


param_update <- function(dp){

  newParams <- ClusterParameterUpdate(dp)

  dp$uniqueParams <- newParams
  dp$params <- newParams[dp$states]
  return(dp)
}
