/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */
#pragma once

#include "SparsifyingRegularizer.h"
#include "OptimizerSparse.h"

using arma::vec;
using arma::uvec;
using arma::uword;
using arma::sp_mat;
using arma::umat;
using Rcpp::List;
using Rcpp::as;
using std::vector;

template <typename matrix, SparseNorm norm>
class SparseRegularizer :
  public SparsifyingRegularizer<matrix> {
  public:

    using SparsifyingRegularizer<matrix>::intercept_ ;
    using SparsifyingRegularizer<matrix>::gamma_     ;
    using SparsifyingRegularizer<matrix>::data_      ;
    using SparsifyingRegularizer<matrix>::beta_      ;
    using SparsifyingRegularizer<matrix>::grad_      ;
    using SparsifyingRegularizer<matrix>::lambda_factor_ ;
    using SparsifyingRegularizer<matrix>::get_lambda_seq ;
    using SparsifyingRegularizer<matrix>::nzeros_    ;
    using SparsifyingRegularizer<matrix>::debiased_  ;
    using SparsifyingRegularizer<matrix>::active_    ;
    using SparsifyingRegularizer<matrix>::intercept_debiased_ ;
    using SparsifyingRegularizer<matrix>::beta_debiased_ ;
    using typename SparsifyingRegularizer<matrix>::StepResult ;
    using typename SparsifyingRegularizer<matrix>::Diagnostics ;

    SparsePenalty<norm>          penalty_ ;
    ActiveSet<matrix>            set_     ;
    SparseOptimizer<matrix,norm> solver_  ;

    double get_lambda_max() {
      return penalty_.lambda_max(data_.XTy_, lambda_factor_) ;
    }

    SparseRegularizer(RegressionData<matrix>, const List&, const List&) ;

    double get_df() override ;

    // ── SparsifyingRegularizer hooks ─────────────────────────────────────────

    ActiveSet<matrix>& current_set() override { return set_ ; }

    StepResult run_solver(double lambda) override {
      uword status = solver_.working_set(beta_, grad_, lambda, lambda_factor_, gamma_, data_, set_) ;
      return { status, solver_.gap_, solver_.iter_ } ;
    }

    Diagnostics solver_diagnostics() const override {
      return { solver_.inner_iter_, solver_.J_vec_, solver_.D_vec_ } ;
    }

} ;

template <typename matrix, SparseNorm norm>
SparseRegularizer<matrix,norm>::SparseRegularizer(
  RegressionData<matrix> data, const List& regParam, const List& control) :
  SparsifyingRegularizer<matrix>::SparsifyingRegularizer(std::move(data), regParam) {

    // Scale the structuring matrix according to the amount of l2 penalty
    data_.scale_struct(gamma_) ;

    // Initialize the active set, beta_ and the gradient
    vec beta0 = control["beta0"] ;
    uvec A0 = find(beta0) ;
    grad_ = - data_.XTy_ ;
    if (A0.is_empty()) {
      set_  = ActiveSet(data_, as<bool>(control["factmat"])) ;
    } else {
      set_  = ActiveSet(data_, A0, as<bool>(control["factmat"])) ;
      beta_ = beta0(A0) ;
      grad_ += set_.XTXA_ * beta_ ;
    }

    // Set the penalty
    penalty_ = SparsePenalty<norm>(as<double>(regParam["eta"])) ;
    get_lambda_seq(penalty_.lambda_max(data_.XTy_, lambda_factor_), regParam) ;

    // Set up the optimizer
    solver_ = SparseOptimizer<matrix,norm>(penalty_, control) ;
  }

template <typename matrix, SparseNorm norm>
double SparseRegularizer<matrix,norm>::get_df() {

  double df = set_.size() + data_.centered_ ;
  if (gamma_ > 0) {
    set_.inverse_Gram() ;
    uword k = set_.size() ;
    // Build local position lookup: pos(i) = index of variable i in set_.A_
    uvec pos(data_.p_) ;
    for (uword i = 0; i < k; i++) pos(set_.A_(i)) = i ;
    // Iterate over non-zeros of S_ only — O(p + nnz(S)) vs O(k² log nnz)
    mat SAA(k, k, arma::fill::zeros) ;
    for (auto it = data_.S_.begin(); it != data_.S_.end(); ++it) {
      uword r = it.row(), c = it.col() ;
      if (set_.is_in_(r) && set_.is_in_(c)) SAA(pos(r), pos(c)) = *it ;
    }
    df -= accu(SAA % set_.XATXAinv_) ;
  }

  return df ;
}
