/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/GroupSparseRegularizer.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List group_enet_l1l2_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<mat,GroupSparseNorm::L1L2> grpenet(std::move(data), group, regParam, control) ;
  return grpenet.to_list(grpenet.solution_path(control)) ;
}

// [[Rcpp::export]]
List group_enet_l1l2_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<sp_mat,GroupSparseNorm::L1L2> grpenet(std::move(data), group, regParam, control) ;
  return grpenet.to_list(grpenet.solution_path(control)) ;
}

// [[Rcpp::export]]
List group_enet_l1linf_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<mat,GroupSparseNorm::L1LINF> grpenet(std::move(data), group, regParam, control) ;
  return grpenet.to_list(grpenet.solution_path(control)) ;
}

// [[Rcpp::export]]
List group_enet_l1linf_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<sp_mat,GroupSparseNorm::L1LINF> grpenet(std::move(data), group, regParam, control) ;
  return grpenet.to_list(grpenet.solution_path(control)) ;
}

// [[Rcpp::export]]
List group_enet_coop_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<mat,GroupSparseNorm::COOP> coop(std::move(data), group, regParam, control) ;
  return coop.to_list(coop.solution_path(control)) ;
}

// [[Rcpp::export]]
List group_enet_coop_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const arma::uvec  &group,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  GroupSparseRegularizer<sp_mat,GroupSparseNorm::COOP> coop(std::move(data), group, regParam, control) ;
  return coop.to_list(coop.solution_path(control)) ;
}
