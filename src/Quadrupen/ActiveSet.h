/*
 * Author: Julien CHIQUET
 * MIA PS
 */
#pragma once

using arma::vec;
using arma::mat;
using arma::uvec;
using arma::uword;
using arma::colvec;
using arma::zeros;
using arma::eye;

#include "RegressionData.h"

template <typename matrix>
class ActiveSet {

public:

  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec A_           ; // set of currently activated variables
  uvec is_in_       ; // indicator of active variables (0/1)
  mat XATXA_, XTXA_ ; // matrices of currently activated variables
  mat XATXAinv_     ;
  bool use_chol_    ; // Maintain a Cholesky factorization along the active set algorithm
  mat R_            ; // Cholesky decomposition of XATXA

  ActiveSet() {} ;
  ActiveSet(const RegressionData<matrix> &data, const bool use_chol=true) ;
  ActiveSet(const RegressionData<matrix> &data, const uvec&, const bool use_chol) ;

  // ── Active set handling ──────────────────────────────────────────────────────────
  void add_var(uword, const RegressionData<matrix> &) ; // add a single variable
  void add_vars(uvec, const RegressionData<matrix> &) ; // add a list of variables
  void del_var(uword, vec&) ; // remove the variable activated in position ind_var_out
  void del_vars(uvec, vec&) ; // remove a set of non contiguous variables
  void reset() ; // empty the active set
  const uword size() const { return A_.n_elem ; }

  // ── Update/Downdate the Cholesky factorisation ────────────────────────────────────
  void update_Cholesky() ; // Insert the last activated variable
  void update_Cholesky_block(uword n_new) ; // Insert the last n_new activated variables
  void downdate_Cholesky(uword j) ; // Remove the specified variables

  // ── Inverse the currently active Gram matrix (XATXAinv_) ──────────────────────────
  void inverse_Gram() ; // When whole inverse is needed (df computation when gamma > 0)
  vec solve_Gram(const vec& b) const ; // Without the full inverse — O(k²) vs O(k³)

};

// ── Constructors ────────────────────────────────────────────────────────────────────
template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
}

template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const uvec& A0, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
  add_vars(A0, data)    ;
}

template <typename matrix>
void ActiveSet<matrix>::reset() {
  A_.reset()      ;
  is_in_.zeros()  ;
  XATXA_.reset()  ;
  XTXA_.reset()   ;
  R_.reset()      ;
}

template <typename matrix>
void ActiveSet<matrix>::add_var(uword var_in, const RegressionData<matrix>& data) {
  uword k = size() ;
  A_.resize(k + 1) ;
  A_(k) = var_in   ;
  is_in_[var_in] = 1 ;

  vec wcol = data.weights_ % vec(data.X_.col(var_in)) ;
  vec new_col = data.X_.t() * wcol -
    data.n_w_ * data.X_bar_ * arma::as_scalar(data.X_bar_[var_in]) + data.S_.col(var_in) ;

  // Single allocation for XTXA_: copy old columns then set new one
  mat new_XTXA_(data.p_, k + 1, arma::fill::none) ;
  if (k > 0) new_XTXA_.cols(0, k - 1) = XTXA_ ;
  new_XTXA_.col(k) = new_col ;
  XTXA_ = std::move(new_XTXA_) ;

  // Single allocation for XATXA_: fill four blocks directly
  // [ XATXA_old | cross        ]
  // [ cross.t() | new_cols.rows(vars) ]
  mat new_XATXA_(k + 1, k + 1, arma::fill::none) ;
  if (k > 0) {
    vec cross = new_col.elem(A_.head(k)) ; // cross-products with previously active variables
    new_XATXA_.submat(0, 0, k-1, k-1) = XATXA_ ;
    new_XATXA_.col(k).head(k)  = cross ;
    new_XATXA_.row(k).head(k)  = cross.t() ;
  }
  new_XATXA_(k, k) = new_col(var_in) ;
  XATXA_ = std::move(new_XATXA_) ;

  if (use_chol_) update_Cholesky() ;
}

template <typename matrix>
void ActiveSet<matrix>::add_vars(uvec vars, const RegressionData<matrix>& data) {
  uword n_new   = vars.n_elem ;
  uword p_old   = size() ;
  uword p_total = p_old + n_new ;

  for (uword v : vars) is_in_[v] = 1 ;
  A_.resize(p_total) ;
  A_.tail(n_new) = vars ;

  mat WXvars(data.X_.cols(vars)) ;
  WXvars.each_col() %= data.weights_ ;
  mat new_cols = data.X_.t() * WXvars -
    data.n_w_ * data.X_bar_ * data.X_bar_.rows(vars).t() +
    data.S_.cols(vars) ;

  // Single allocation for XTXA_
  mat new_XTXA_(data.p_, p_total, arma::fill::none) ;
  if (p_old > 0) new_XTXA_.cols(0, p_old - 1) = XTXA_ ;
  new_XTXA_.cols(p_old, p_total - 1) = new_cols ;
  XTXA_ = std::move(new_XTXA_) ;

  // Single allocation for XATXA_
  mat new_XATXA_(p_total, p_total, arma::fill::none) ;
  if (p_old > 0) {
    mat cross = new_cols.rows(A_.head(p_old)) ; // p_old x n_new cross-products
    new_XATXA_.submat(0,     0,     p_old-1,   p_old-1)   = XATXA_ ;
    new_XATXA_.submat(0,     p_old, p_old-1,   p_total-1) = cross ;
    new_XATXA_.submat(p_old, 0,     p_total-1, p_old-1)   = cross.t() ;
  }
  new_XATXA_.submat(p_old, p_old, p_total-1, p_total-1) = new_cols.rows(vars) ;
  XATXA_ = std::move(new_XATXA_) ;

  if (use_chol_) update_Cholesky_block(n_new) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_var(uword ivar_out, vec& beta) {
  is_in_[A_[ivar_out]] = 0  ;
  A_.shed_row(ivar_out)     ;
  XTXA_.shed_col(ivar_out)  ;
  XATXA_.shed_col(ivar_out) ;
  XATXA_.shed_row(ivar_out) ;
  beta.shed_row(ivar_out)   ;

  if (use_chol_) downdate_Cholesky(ivar_out) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_vars(uvec ivars, vec& beta) {
  ivars = sort(ivars, "descend");
  for (uword i=0 ; i <ivars.n_elem ; i++) {
    del_var(ivars[i], beta) ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky() {
  uword p = XATXA_.n_cols ;

  if (p == 1) {
    R_ = sqrt(XATXA_) ;
  } else {
    // Solve R_old^T * rp[0..p-2] = XATXA_[0..p-2, p-1]
    colvec rp(p, arma::fill::zeros) ;
    rp.head(p-1) = solve(trimatu(R_).t(),
                         XATXA_.col(p-1).head(p-1),
                         arma::solve_opts::fast) ;
    rp(p-1) = std::sqrt(XATXA_(p-1, p-1) - dot(rp.head(p-1), rp.head(p-1))) ;

    // Extend R_ from (p-1)x(p-1) to pxp
    // [ R_old | R_new_cols     ]
    // [ 0     | R_bottom_right ]
    mat new_R_(p, p, arma::fill::zeros) ; // lower-triangular part stays zero
    new_R_.submat(0, 0, p-2, p-2) = R_ ;
    new_R_.col(p-1) = rp ;
    R_ = std::move(new_R_) ;
  }
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky_block(uword n_new) {
  uword p_total = XATXA_.n_cols ;
  uword p_old   = p_total - n_new ;

  if (p_old == 0) {
    R_ = chol(XATXA_) ;
  } else {
    // Solve R_old^T * R_new_cols = XATXA_[0..p_old-1, p_old..p_total-1]
    mat R_new_cols = solve(trimatu(R_).t(),
                           XATXA_.submat(0, p_old, p_old-1, p_total-1),
                           arma::solve_opts::fast) ;

    // Schur complement for the new diagonal block
    mat R_bottom_right = chol(XATXA_.submat(p_old, p_old, p_total-1, p_total-1) -
                              R_new_cols.t() * R_new_cols) ;

    // Extend R_ from p_old×p_old to p_total×p_total
    // [ R_old | R_new_cols     ]
    // [ 0     | R_bottom_right ]
    mat new_R_(p_total, p_total, arma::fill::zeros) ; // lower-triangular part stays zero
    new_R_.submat(0,     0,     p_old-1,   p_old-1)   = R_ ;
    new_R_.submat(0,     p_old, p_old-1,   p_total-1) = R_new_cols ;
    new_R_.submat(p_old, p_old, p_total-1, p_total-1) = R_bottom_right ;
    R_ = std::move(new_R_) ;
  }
}
template <typename matrix>
void ActiveSet<matrix>::downdate_Cholesky(uword j) {

  vec x = zeros<vec>(2);
  mat G = zeros<mat>(2,2);

  R_.shed_col(j);
  int p = R_.n_cols;
  double r;
  for (int k=j; k<p; k++) {
    x = R_.submat(k,k,k+1,k);

    if (x[1] != 0) {
      r = std::hypot(x(0), x(1));
      G = {{x(0), x(1)}, {-x(1), x(0)}};
      G = G / r;
      x(0) = r; x(1) = 0;
    } else {
      G = eye(2,2);
    }
    R_.submat(k,k,k+1,k) = x;
    if (k < p-1) {
      R_.submat(k,k+1,k+1,p-1) = G * R_.submat(k,k+1,k+1,p-1);
    }
  }
  R_.shed_row(p);
}

template <typename matrix>
void ActiveSet<matrix>::inverse_Gram() {
  if (use_chol_) {
    XATXAinv_ = solve(trimatu(R_), solve(trimatl(R_.t()), eye<mat>(R_.n_cols, R_.n_cols)));
  } else {
    XATXAinv_ = inv_sympd(XATXA_, arma::inv_opts::allow_approx);
  }
}

template <typename matrix>
vec ActiveSet<matrix>::solve_Gram(const vec& b) const {
  if (use_chol_) {
    return solve(trimatu(R_),
                 solve(trimatl(R_.t()), b, arma::solve_opts::fast),
                 arma::solve_opts::fast) ;
  } else {
    return solve(XATXA_, b, arma::solve_opts::fast) ;
  }
}
