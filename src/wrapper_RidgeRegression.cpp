/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/RidgeRegression.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List ridge_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &controlFit
    ) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(controlFit["normalize"])) ;
  RidgeRegression ridge(std::move(data), regParam) ;
  return ridge.to_list(ridge.solution_path(trimatu(as<mat>(dataModel["C_inv"])))) ;
}
