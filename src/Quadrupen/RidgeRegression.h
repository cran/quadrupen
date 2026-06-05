/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "Regularizer.h"
#include "PenaltyDense.h"

using arma::mat;
using Rcpp::List;

class RidgeRegression : public Regularizer<mat>{
public:
  RidgeRegression(RegressionData<mat>, const List&);

  double get_lambda_max() {return(penalty_.dual_norm(data_.XTy_, lambda_factor_));}
  
  List solution_path(const mat&);
  List to_list(const List& monitoring) ;
  
  DensePenalty<DenseNorm::L2> penalty_ ; // main penalty object 
  
};

