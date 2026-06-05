/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */
#pragma once

#include "Regularizer.h"
#include "PenaltyDense.h"
#include "OptimizerLINF.h"

using arma::mat;
using arma::sp_mat;
using arma::uvec;
using arma::umat;
using arma::uword;
using arma::vec;
using Rcpp::List;
using std::vector;

class BoundedRegression : public Regularizer<mat> {
public:
  
  // Specific to Bounded regression
  DensePenalty<DenseNorm::LINF> penalty_ ; // main penalty object 
  OptimizerLINF solver_ ; // Solvers for LINF penalty
  uvec unbounded_   ; // Active variables (away from the boundary)
  vector<uvec> bounded_ ; // variables reaching the boundary (for all lambda values)
  
  BoundedRegression(RegressionData<mat>, const List&, const List&);
  
  double get_lambda_max() {
    return(penalty_.dual_norm(data_.XTy_, lambda_factor_));
  }
  
  const sp_mat unbounded_var() {
    umat locs = build_sp_locations(bounded_) ;
    return sp_mat(locs, arma::ones<vec>(locs.n_cols),
                  data_.p_, bounded_.size(), true, false) ;
  }
  
  List solution_path(const List&);
  List to_list(const List& monitoring) ;
  
  // Compute degrees of freedom for the current estimate
  double get_df() ; 
  
};


