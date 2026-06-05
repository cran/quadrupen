/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupLava.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_lava_l1l2_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupLava<mat,GroupSparseNorm::L1L2> group_lava(std::move(data), group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;
  return group_lava.to_list(results) ;
}

// [[Rcpp::export]]
List group_lava_l1linf_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupLava<mat,GroupSparseNorm::L1LINF> group_lava(std::move(data), group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;
  return group_lava.to_list(results) ;
}

// [[Rcpp::export]]
List group_lava_coop_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupLava<mat,GroupSparseNorm::COOP> group_lava(std::move(data), group, regParam, control) ;
  List results = group_lava.solution_path(control) ;
  group_lava.post_treatment() ;
  return group_lava.to_list(results) ;
}
