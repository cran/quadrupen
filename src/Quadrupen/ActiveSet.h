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
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

// The k x k matrices of the active set (Gram matrix, Cholesky factor) grow and shrink by a few
// rows and columns at a time. To avoid a reallocation and a full copy at each change, they keep
// spare capacity: Armadillo reuses the allocated memory as long as the number of elements does
// not exceed n_alloc, and the columns are moved in place to the new layout.
namespace square_inplace {

inline void check_same_memory(const mat& M, const double* mem) {
  if (M.memptr() != mem) Rcpp::stop("internal error: unexpected reallocation of an active set matrix") ;
}

// Grow M from k x k to (k+m) x (k+m), keeping its top-left block; new entries are not initialized
inline void grow(mat& M, uword m) {
  const uword k = M.n_rows, kn = k + m ;
  if (kn * kn > M.n_alloc) {
    // not enough capacity: reallocate with ~56% spare elements (side x 1.25)
    const uword cap = std::max(kn, (uword) std::ceil(1.25 * kn)) ;
    mat bigger(cap, cap, arma::fill::none) ;
    bigger.set_size(kn, kn) ; // keeps the cap x cap allocation
    for (uword j = 0; j < k; ++j) std::copy_n(M.colptr(j), k, bigger.colptr(j)) ;
    M = std::move(bigger) ;
    return ;
  }
  // move the columns from leading dimension k to kn, last column first (destinations are after sources)
  double* mem = M.memptr() ;
  for (uword j = k; j-- > 1; ) std::memmove(mem + j * kn, mem + j * k, k * sizeof(double)) ;
  M.set_size(kn, kn) ;
  check_same_memory(M, mem) ;
}

// Remove row and column i of the square matrix M
inline void remove(mat& M, uword i) {
  const uword k = M.n_rows, kn = k - 1 ;
  if (kn * kn <= arma::arma_config::mat_prealloc) { M.shed_col(i) ; M.shed_row(i) ; return ; }
  // compact the kept entries in column-major order (destinations are never after sources)
  double* mem = M.memptr() ;
  uword dst = 0 ;
  for (uword j = 0; j < k; ++j) {
    if (j == i) continue ;
    const double* col = mem + j * k ;
    for (uword r = 0; r < k; ++r) if (r != i) mem[dst++] = col[r] ;
  }
  M.set_size(kn, kn) ;
  check_same_memory(M, mem) ;
}

} // namespace square_inplace

template <typename matrix>
class ActiveSet {

public:

  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec A_           ; // set of currently activated variables
  uvec is_in_       ; // indicator of active variables (0/1)
  mat XATXA_        ; // Gram matrix of currently activated variables
  mat XATXAinv_     ;
  bool use_chol_    ; // Maintain a Cholesky factorization along the active set algorithm
  mat R_            ; // Cholesky decomposition of XATXA (upper triangular, XATXA = R'R)

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

  // ── Products with X'X_A (p x k), stored without copy ──────────────────────────────
  vec XTXA_times(const vec& v) const ;

  // ── Update/Downdate the Cholesky factorisation ────────────────────────────────────
  void update_Cholesky() ; // Insert the last activated variable
  void update_Cholesky_block(uword n_new) ; // Insert the last n_new activated variables
  void downdate_Cholesky(uword j) ; // Remove the specified variables

  // ── Inverse the currently active Gram matrix (XATXAinv_) ──────────────────────────
  void inverse_Gram() ; // When whole inverse is needed (df computation when gamma > 0)
  vec solve_Gram(const vec& b) const ; // Without the full inverse — O(k²) vs O(k³)

private:

  // The first size() columns of XTXA_buf_ hold X'X_A. Extra columns are spare
  // capacity, so that adding/removing a variable does not reallocate a p x k matrix.
  mat XTXA_buf_ ;
  void reserve_XTXA(uword k_old, uword k_new) ;
  void compact_XTXA(const uvec& positions, uword k_old) ;

  // Extend R_ from (p-1) x (p-1) to p x p with the p-th column of XATXA_
  void extend_Cholesky(uword p) ;

  // In-place triangular solve with R_: R' X = B (trans = 'T') or R X = B (trans = 'N')
  // Returns false when R_ is singular or the solution is not finite
  bool solve_R(mat& B, char trans) const ;
};

// ── Constructors ────────────────────────────────────────────────────────────────────
template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
  XTXA_buf_.set_size(data.p_, 0) ;
}

template <typename matrix>
ActiveSet<matrix>::ActiveSet(const RegressionData<matrix>& data, const uvec& A0, const bool use_chol) :
  use_chol_(use_chol) {
  is_in_.zeros(data.p_) ;
  XTXA_buf_.set_size(data.p_, 0) ;
  add_vars(A0, data)    ;
}

template <typename matrix>
void ActiveSet<matrix>::reset() {
  A_.reset()      ;
  is_in_.zeros()  ;
  XATXA_.reset()  ;
  XTXA_buf_.set_size(XTXA_buf_.n_rows, 0) ;
  R_.reset()      ;
}

// ── X'X_A storage ───────────────────────────────────────────────────────────────────
template <typename matrix>
vec ActiveSet<matrix>::XTXA_times(const vec& v) const {
  if (size() == 0) return zeros<vec>(XTXA_buf_.n_rows) ;
  // read-only alias on the first size() columns of the buffer (no copy)
  const mat XTXA(const_cast<double*>(XTXA_buf_.memptr()), XTXA_buf_.n_rows, size(), false, true) ;
  return XTXA * v ;
}

template <typename matrix>
void ActiveSet<matrix>::reserve_XTXA(uword k_old, uword k_new) {
  if (k_new <= XTXA_buf_.n_cols) return ;
  uword p = XTXA_buf_.n_rows ;
  // geometric growth, bounded by p (the active set cannot exceed p variables)
  uword capacity = std::max(k_new, std::min(p, std::max<uword>(16, 2 * XTXA_buf_.n_cols))) ;
  mat new_buf(p, capacity, arma::fill::none) ;
  if (k_old > 0) std::copy_n(XTXA_buf_.memptr(), p * k_old, new_buf.memptr()) ;
  XTXA_buf_ = std::move(new_buf) ;
}

template <typename matrix>
void ActiveSet<matrix>::compact_XTXA(const uvec& positions, uword k_old) {
  // Remove the columns at 'positions' among the first k_old ones, preserving order
  std::vector<bool> removed(k_old, false) ;
  for (uword i : positions) removed[i] = true ;
  uword p = XTXA_buf_.n_rows, dst = 0 ;
  for (uword j = 0; j < k_old; ++j) {
    if (removed[j]) continue ;
    if (dst != j) std::copy_n(XTXA_buf_.colptr(j), p, XTXA_buf_.colptr(dst)) ;
    ++dst ;
  }
}

// ── Active set handling ─────────────────────────────────────────────────────────────
template <typename matrix>
void ActiveSet<matrix>::add_var(uword var_in, const RegressionData<matrix>& data) {
  uword k = size() ;
  A_.resize(k + 1) ;
  A_(k) = var_in   ;
  is_in_[var_in] = 1 ;

  vec wcol = data.weights_ % vec(data.X_.col(var_in)) ;
  vec new_col = data.X_.t() * wcol -
    data.n_w_ * data.X_bar_ * arma::as_scalar(data.X_bar_[var_in]) + data.S_.col(var_in) ;

  reserve_XTXA(k, k + 1) ;
  XTXA_buf_.col(k) = new_col ;

  // Grow XATXA_ in place and fill the new row and column
  // [ XATXA_old | cross       ]
  // [ cross.t() | new_col(j)  ]
  square_inplace::grow(XATXA_, 1) ;
  if (k > 0) {
    const vec cross = new_col.elem(A_.head(k)) ; // cross-products with previously active variables
    XATXA_.col(k).head(k) = cross ;
    XATXA_.row(k).head(k) = cross.t() ;
  }
  XATXA_(k, k) = new_col(var_in) ;

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

  reserve_XTXA(p_old, p_total) ;
  XTXA_buf_.cols(p_old, p_total - 1) = new_cols ;

  // Grow XATXA_ in place and fill the new blocks
  square_inplace::grow(XATXA_, n_new) ;
  if (p_old > 0) {
    const mat cross = new_cols.rows(A_.head(p_old)) ; // p_old x n_new cross-products
    XATXA_.submat(0,     p_old, p_old-1,   p_total-1) = cross ;
    XATXA_.submat(p_old, 0,     p_total-1, p_old-1)   = cross.t() ;
  }
  XATXA_.submat(p_old, p_old, p_total-1, p_total-1) = new_cols.rows(vars) ;

  if (use_chol_) update_Cholesky_block(n_new) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_var(uword ivar_out, vec& beta) {
  compact_XTXA(uvec{ivar_out}, size()) ;
  is_in_[A_[ivar_out]] = 0  ;
  A_.shed_row(ivar_out)     ;
  square_inplace::remove(XATXA_, ivar_out) ;
  beta.shed_row(ivar_out)   ;

  if (use_chol_) downdate_Cholesky(ivar_out) ;
}

template <typename matrix>
void ActiveSet<matrix>::del_vars(uvec ivars, vec& beta) {
  if (ivars.is_empty()) return ;
  compact_XTXA(ivars, size()) ; // single pass over X'X_A for all removed variables
  ivars = sort(ivars, "descend");
  for (uword i=0 ; i <ivars.n_elem ; i++) {
    uword ivar_out = ivars[i] ;
    is_in_[A_[ivar_out]] = 0  ;
    A_.shed_row(ivar_out)     ;
    square_inplace::remove(XATXA_, ivar_out) ;
    beta.shed_row(ivar_out)   ;
    if (use_chol_) downdate_Cholesky(ivar_out) ;
  }
}

// ── Cholesky factorisation ──────────────────────────────────────────────────────────
template <typename matrix>
bool ActiveSet<matrix>::solve_R(mat& B, char trans) const {
  // LAPACK dtrtrs directly on R_ memory: no transpose, no copy, no rcond estimation
  if (R_.n_rows == 0 || B.n_elem == 0) return true ;
  char uplo = 'U', diag = 'N' ;
  arma::blas_int n = R_.n_rows, nrhs = B.n_cols, info = 0 ;
  arma::lapack::trtrs<double>(&uplo, &trans, &diag, &n, &nrhs, R_.memptr(), &n, B.memptr(), &n, &info) ;
  return (info == 0) && B.is_finite() ;
}

template <typename matrix>
void ActiveSet<matrix>::extend_Cholesky(uword p) {
  if (p == 1) {
    R_ = sqrt(XATXA_.submat(0, 0, 0, 0)) ;
    return ;
  }
  // Solve R_old^T * rp = XATXA_[0..p-2, p-1]
  vec rp = XATXA_.col(p-1).head(p-1) ;
  solve_R(rp, 'T') ;

  // Extend R_ from (p-1)x(p-1) to pxp in place
  // [ R_old | rp             ]
  // [ 0     | R_bottom_right ]
  // When XATXA_ is (numerically) singular the new pivot is not finite: solve_Gram then
  // falls back to a direct solve.
  square_inplace::grow(R_, 1) ;
  R_.row(p-1).head(p-1).zeros() ;
  R_.col(p-1).head(p-1) = rp ;
  R_(p-1, p-1) = std::sqrt(XATXA_(p-1, p-1) - dot(rp, rp)) ;
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky() {
  extend_Cholesky(XATXA_.n_cols) ;
}

template <typename matrix>
void ActiveSet<matrix>::update_Cholesky_block(uword n_new) {
  uword p_total = XATXA_.n_cols ;
  uword p_old   = p_total - n_new ;

  if (p_old == 0) {
    if (arma::chol(R_, XATXA_)) return ;
  } else {
    // Solve R_old^T * R_new_cols = XATXA_[0..p_old-1, p_old..p_total-1]
    mat R_new_cols = XATXA_.submat(0, p_old, p_old-1, p_total-1) ;
    solve_R(R_new_cols, 'T') ;

    // Schur complement for the new diagonal block
    mat R_bottom_right ;
    if (arma::chol(R_bottom_right, XATXA_.submat(p_old, p_old, p_total-1, p_total-1) -
                                   R_new_cols.t() * R_new_cols)) {
      // Extend R_ from p_old×p_old to p_total×p_total in place
      // [ R_old | R_new_cols     ]
      // [ 0     | R_bottom_right ]
      square_inplace::grow(R_, n_new) ;
      R_.submat(p_old, 0, p_total-1, p_old-1).zeros() ;
      R_.submat(0,     p_old, p_old-1,   p_total-1) = R_new_cols ;
      R_.submat(p_old, p_old, p_total-1, p_total-1) = R_bottom_right ;
      return ;
    }
  }
  // The new block is not positive definite (e.g. more active variables than the rank of X):
  // insert the variables one at a time, as for single additions
  for (uword p = p_old + 1; p <= p_total; ++p) extend_Cholesky(p) ;
}

template <typename matrix>
void ActiveSet<matrix>::downdate_Cholesky(uword j) {
  // Remove column j, then restore the upper triangular form with Givens rotations
  // applied in place on rows (k, k+1), and drop the last row. All steps work in the
  // memory of R_ with its current leading dimension n.
  const uword n = R_.n_rows, p = n - 1 ;
  if (p * p <= arma::arma_config::mat_prealloc) {
    R_.shed_col(j) ;
    for (uword k = j; k < p; ++k) {
      const double a = R_(k, k), b = R_(k+1, k) ;
      if (b == 0.0) continue ;
      const double r = std::hypot(a, b), c = a / r, s = b / r ;
      R_(k, k) = r ; R_(k+1, k) = 0.0 ;
      for (uword l = k + 1; l < p; ++l) {
        const double x = R_(k, l), y = R_(k+1, l) ;
        R_(k, l) = c * x + s * y ; R_(k+1, l) = c * y - s * x ;
      }
    }
    R_.shed_row(p) ;
    return ;
  }

  double* mem = R_.memptr() ;
  std::memmove(mem + j * n, mem + (j + 1) * n, (p - j) * n * sizeof(double)) ; // drop column j
  for (uword k = j; k < p; ++k) {
    double* ck = mem + k * n ;
    const double a = ck[k], b = ck[k+1] ;
    if (b == 0.0) continue ;
    const double r = std::hypot(a, b), c = a / r, s = b / r ;
    ck[k] = r ; ck[k+1] = 0.0 ;
    for (uword l = k + 1; l < p; ++l) {
      double* cl = mem + l * n ;
      const double x = cl[k], y = cl[k+1] ;
      cl[k]   = c * x + s * y ;
      cl[k+1] = c * y - s * x ;
    }
  }
  for (uword l = 1; l < p; ++l) std::memmove(mem + l * p, mem + l * n, p * sizeof(double)) ; // drop the last row
  R_.set_size(p, p) ;
  square_inplace::check_same_memory(R_, mem) ;
}

template <typename matrix>
void ActiveSet<matrix>::inverse_Gram() {
  if (use_chol_) {
    XATXAinv_ = eye<mat>(R_.n_cols, R_.n_cols) ;
    if (solve_R(XATXAinv_, 'T') && solve_R(XATXAinv_, 'N')) return ;
  }
  // no factorization, or degenerate one: approximate inverse
  XATXAinv_ = inv_sympd(XATXA_, arma::inv_opts::allow_approx);
}

template <typename matrix>
vec ActiveSet<matrix>::solve_Gram(const vec& b) const {
  if (use_chol_) {
    vec x = b ;
    if (solve_R(x, 'T') && solve_R(x, 'N')) return x ;
  }
  // no factorization, or degenerate one: fall back to a direct solve, and to a least squares
  // solution when XATXA_ is singular
  vec x ;
  if (arma::solve(x, XATXA_, b, arma::solve_opts::fast)) return x ;
  return arma::solve(XATXA_, b, arma::solve_opts::force_approx) ;
}
