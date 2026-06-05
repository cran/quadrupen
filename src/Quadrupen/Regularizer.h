/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#pragma once

#include "RegressionData.h"

using arma::vec;
using arma::uvec;
using arma::umat;
using arma::uword;
using arma::conv_to;
using arma::logspace;
using Rcpp::List;
using Rcpp::as;
using std::vector;

template <typename matrix>
class Regularizer {
public:

  Regularizer() {} ;
  Regularizer(RegressionData<matrix>, const List&);

  RegressionData<matrix> data_ ; // data structure
  double gamma_ = 0.0          ; // overall amount of minor penalty (not leading the path)
  vector<double> lambdas_      ; // vector of parameters tuning the main penalty
  vec lambda_factor_           ; // weights for the main penalty
  vector<double> intercept_    ; // vector of intercept term
  matrix coef_                 ; // matrix of coefficients
  vec beta_                    ; // vector of current parameters (for fix lambda value)
  vec grad_                    ; // vector of current gradient (smooth part)
  vector<double> df_           ; // degrees of freedom along the path

  void get_lambda_seq(double, const List&);

  // Build (row, col) location matrix for sp_mat construction from a sequence
  // of active index sets. Row 0 = variable indices, Row 1 = lambda-step index.
  static umat build_sp_locations(const vector<uvec>& groups) {
    uword total_nnz = 0;
    for (const auto& g : groups) total_nnz += g.n_elem;
    umat locs(2, total_nnz);
    uword col = 0, pos = 0;
    for (const auto& g : groups) {
      for (uword r : g) { locs(0, pos) = r; locs(1, pos) = col; ++pos; }
      ++col;
    }
    return locs;
  }
};

template <typename matrix>
Regularizer<matrix>::Regularizer(
  RegressionData<matrix> data, const List& regParam) :
  data_ (std::move(data)),
  gamma_(as<double>(regParam["gamma"])),
  lambda_factor_(as<vec>(regParam["lambda_factor"]))
{}

template <typename matrix>
void Regularizer<matrix>::get_lambda_seq(double lambda_max, const List& regParam) {
  if (regParam["lambda"] != R_NilValue) {
    lambdas_  = as<vector<double>>(regParam["lambda"]) ;
  } else {
    lambdas_ = conv_to<vector<double>>::from(
      logspace(
        log10(lambda_max),
        log10(as<double>(regParam["min_ratio"])*lambda_max),
        as<uword>(regParam["n_lambda"])
      )
    );
  }
}
