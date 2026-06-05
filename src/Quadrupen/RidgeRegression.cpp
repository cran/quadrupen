/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "RidgeRegression.h"

using namespace Rcpp;
using namespace arma;

RidgeRegression::RidgeRegression(
  RegressionData<mat> data, const List& regParam) :
  Regularizer::Regularizer(std::move(data), regParam)
{
  penalty_ = DensePenalty<DenseNorm::L2>() ;
  get_lambda_seq(get_lambda_max(), regParam) ;
}

List RidgeRegression::to_list(const List& monitoring) {
  return List::create(
    Named("tuning_param") = List::create(
      Named("l2") = lambdas_,
      Named("l1") = 0
    ),
    Named("coef")       = coef_,
    Named("intercept")  = intercept_,
    Named("normx")      = data_.norm_X_,
    Named("df")         = df_,
    Named("monitoring") = monitoring
  ) ;
}

List RidgeRegression::solution_path(const mat& C_inv) {

  // Weighted SVD: diag(sqrtW) * Xc * C_inv
  vec sqrtW = sqrt(data_.weights_) ;
  mat Xwc   = data_.X_ ;
  Xwc.each_row() -= data_.X_bar_.t() ;
  Xwc.each_col() %= sqrtW ;

  vec eta ; mat U, V ;
  svd_econ(U, eta, V, Xwc * C_inv) ;

  mat C_invV = C_inv * V ;
  vec Uty    = U.t() * (sqrtW % (data_.y_ - data_.y_bar_)) ;

  vector<double> timing ; // successive timing for solving for each lambda value
  wall_clock timer ; timer.tic(); // clock
  for(auto lambda : lambdas_) {
    // computing the structured ridge estimate
    beta_ = (C_invV * diagmat(eta/(square(eta) + lambda)) * Uty) ;
    coef_ = join_rows(coef_, beta_ / data_.norm_X_) ;
    // estimating the intercept term
    intercept_.push_back(data_.y_bar_ - dot(beta_, data_.X_bar_));  
    // computing the estimated degrees of freedom
    df_.push_back(sum(square(eta)/(square(eta) + lambda)) + data_.centered_);
    
    timing.push_back(timer.toc()) ;
  }

  return(
    List::create(
      Named("max_grd")        = R_NaReal,
      Named("convergence")    = 0,
      Named("pensteps_timer") = timing 
    )
  );
  
}
