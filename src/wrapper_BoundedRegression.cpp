/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/BoundedRegression.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List bounded_regression_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
    ) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  BoundedRegression bounded(std::move(data), regParam, control) ;
  return bounded.to_list(bounded.solution_path(control)) ;
}
