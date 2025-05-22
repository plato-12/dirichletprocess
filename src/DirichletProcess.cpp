// src/DirichletProcess.cpp
#include "../inst/include/DirichletProcessBase.h"
#include "../inst/include/RcppConversions.h"

namespace dp {

DirichletProcess::DirichletProcess() : n(0), alpha(1.0) {
  // Default constructor
}

DirichletProcess::~DirichletProcess() {
  // Destructor
}

Rcpp::List DirichletProcess::toR() const {
  Rcpp::List result;

  // Convert common properties
  result["data"] = Rcpp::wrap(data);
  result["n"] = n;
  result["alpha"] = alpha;
  result["alphaPriorParameters"] = alphaPriorParameters;

  // Add class attributes (to be extended by derived classes)
  result.attr("class") = Rcpp::CharacterVector::create("list", "dirichletprocess");

  return result;
}

DirichletProcess* DirichletProcess::fromR(const Rcpp::List& rObj) {
  // This is a factory method that should be implemented by createDPFromR
  return createDPFromR(rObj);
}

MixingDistribution::MixingDistribution() : conjugate(true) {
  // Default constructor
}

MixingDistribution::~MixingDistribution() {
  // Destructor
}

Rcpp::List MixingDistribution::toR() const {
  Rcpp::List result;

  result["distribution"] = distribution;
  result["priorParameters"] = priorParameters;
  result["conjugate"] = conjugate;

  // Use proper Rcpp null checking
  if (mhStepSize.sexp_type() != NILSXP) {
    result["mhStepSize"] = mhStepSize;
  }

  if (hyperPriorParameters.sexp_type() != NILSXP) {
    result["hyperPriorParameters"] = hyperPriorParameters;
  }

  result.attr("class") = Rcpp::CharacterVector::create(
    "list", "MixingDistribution", distribution,
    conjugate ? "conjugate" : "nonconjugate");

  return result;
}

} // namespace dp
