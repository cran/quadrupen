/*
 * Author: Julien CHIQUET, INRAE
 *         julien.chiquet@inrae.fr
 *         Statistique et Génome
 *         MIA Paris-Saclay
 */

#pragma once

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#include <RcppArmadillo.h>

using arma::uword;
using arma::sp_mat;
using arma::vec;
using arma::mat;
using arma::zeros;
using arma::ones;
using Rcpp::Environment;
using Rcpp::as;

// Use template to handle dense or sparse encoding (mat/sp_mat in armadillo)
template <typename matrix>
class RegressionData {

public:

  // DATA VARIABLES FOR REGRESSION PURPOSE
  uword  n_        ; // sample size
  uword  p_        ; // # of features
  matrix X_        ; // matrix of predictors
  vec    y_        ; // vector of response
  sp_mat S_        ; // Structuring matrix
  vec    weights_  ; // observation weights
  double n_w_ = 0.0; // sum of weights (= n_ when weights are uniform)
  bool   centered_ ; // should intercept be considered?
  bool   scaled_   ; // should predictors be standardized?
  mat    XTX_      ; // Gram matrix
  vec    XTy_      ; // responses to predictors vector
  vec    X_bar_    ; // weighted mean of the predictors
  vec    norm_X_   ; // weighted L2 norm of the predictors (centered)
  double y_bar_    ; // weighted mean of the response
  double norm_y_   ; // weighted L2 norm of the response
  
  RegressionData() ;
  RegressionData(const Environment& dataModel, const bool& center, const bool& scale);
  RegressionData(const matrix X,  const vec y, const sp_mat S, const vec weights, const bool center, const bool scale) ;

  void standardize();

  void scale_struct(const double gamma) ;

  void precompute_XTX() ;
    
};


template <typename matrix>
RegressionData<matrix>::RegressionData(
  const Environment& dataModel,
  const bool& center, const bool& scale):
  y_       (as<vec>    (dataModel["y"])) , // response vector
  S_       (as<sp_mat> (dataModel["S"])) , // structuring matrix
  weights_ (as<vec>    (dataModel["wy"])), // observation weights
  centered_(center), scaled_(scale)        // standardization options
{
  X_ = as<matrix>(dataModel["X"]) ;
  n_ = X_.n_rows ; p_ = X_.n_cols ;
  standardize() ;
}

// Constructor from arma objects
template <typename matrix>
RegressionData<matrix>::RegressionData(
  const matrix X, 
  const vec y, 
  const sp_mat S, 
  const vec weights,
  const bool center, const bool scale):
  X_       (X)    , // response vector
  y_       (y) , // response vector
  S_       (S) , // structuring matrix
  weights_ (weights), // observation weights
  centered_(center), scaled_(scale) // standardization options
{ 
  n_ = X_.n_rows ; p_ = X_.n_cols ; 
  standardize() ;
}

template <typename matrix>
void RegressionData<matrix>::scale_struct(const double gamma) {
  S_ *= gamma ; 
};

template <typename matrix>
void RegressionData<matrix>::precompute_XTX() {
  mat WX = X_ ;
  WX.each_col() %= weights_ ;
  XTX_ = X_.t() * WX - n_w_ * X_bar_ * X_bar_.t() + S_ ;
};

template <typename matrix>
void RegressionData<matrix>::standardize() {

  n_w_ = sum(weights_) ;

  if (centered_) {
    X_bar_ = X_.t() * weights_ / n_w_ ;
    y_bar_ = dot(weights_, y_) / n_w_ ;
  } else {
    X_bar_ = zeros(p_) ;
    y_bar_ = 0 ;
  }

  if (scaled_) {
    norm_X_ = sqrt((weights_.t() * (X_ % X_)).t() - n_w_ * square(X_bar_)) ;
    norm_X_.replace(0.0, 1.0) ;
    X_.each_row() /= norm_X_.t() ;
    X_bar_ /= norm_X_ ;
  } else {
    norm_X_ = ones(p_) ;
  }
  norm_y_ = sqrt(dot(weights_, square(y_))) ;

  XTy_ = X_.t() * (weights_ % (y_ - y_bar_)) ;
}

// Specialization for sp_mat: column scaling via non-zero iterator — O(nnz) instead
// of the dense each_col() path which would densify the matrix.
// The generic path uses X_.each_col() /= norm_X_ which is O(n*p) for dense mat.
template <>
inline void RegressionData<sp_mat>::standardize() {

  n_w_ = sum(weights_) ;

  if (centered_) {
    X_bar_ = X_.t() * weights_ / n_w_ ;
    y_bar_ = dot(weights_, y_) / n_w_ ;
  } else {
    X_bar_ = zeros(p_) ;
    y_bar_ = 0 ;
  }

  if (scaled_) {
    // weighted column squared norms via non-zero iterator — O(nnz)
    vec col_sq_wt = zeros(p_) ;
    for (auto it = X_.begin() ; it != X_.end() ; ++it)
      col_sq_wt[it.col()] += weights_[it.row()] * (*it) * (*it) ;
    norm_X_ = sqrt(col_sq_wt - n_w_ * square(X_bar_)) ;
    norm_X_.replace(0.0, 1.0) ;
    for (auto it = X_.begin() ; it != X_.end() ; ++it)
      *it /= norm_X_[it.col()] ;
    X_bar_ /= norm_X_ ;
  } else {
    norm_X_ = ones(p_) ;
  }
  norm_y_ = sqrt(dot(weights_, square(y_))) ;

  XTy_ = X_.t() * (weights_ % (y_ - y_bar_)) ;
}
