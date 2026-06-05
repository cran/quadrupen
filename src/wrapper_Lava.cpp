/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/Lava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List lava_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  Lava<mat,SparseNorm::L1> lava(std::move(data), regParam, control) ;
  List results = lava.solution_path(control) ;
  lava.post_treatment() ;
  return lava.to_list(results) ;
}
