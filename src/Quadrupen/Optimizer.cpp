/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Optimizer.h"

using namespace Rcpp;
using namespace arma;

Optimizer::Optimizer(const List& control) :
  accuracy_(control["threshold"]),
verbosity_(control["verbose"]),
maxiter_(control["maxiter"]),
maxfeat_(control["maxfeat"]),
monitoring_(control["monitor"]) {
  
  if (as<std::string>(control["method"]) == "FISTA") algorithm_ = SolverType::FISTA;
  if (as<std::string>(control["method"]) == "QUADRA") algorithm_ = SolverType::QUADRA;
  if (as<std::string>(control["method"]) == "PGD") algorithm_ = SolverType::PGD;
  if (control.containsElementNamed("maxadd")) max_add_ = std::max<uword>(1, as<uword>(control["maxadd"]));
}

uvec Optimizer::select_violators(
  const vec& optimality,
  const uvec& is_in,
  const double& tol,
  const uword& max_add) const {

  uvec candidates = find(optimality > tol && is_in == 0) ;
  if (candidates.n_elem > max_add) {
    uvec order = sort_index(optimality(candidates), "descend") ;
    candidates = candidates(order.head(max_add)) ;
  } else if (candidates.n_elem > 1) {
    candidates = candidates(sort_index(optimality(candidates), "descend")) ;
  }
  return candidates ;
}

double Optimizer::estimate_lipschitz(
  const mat& XTX,
  uword max_it,
  double tol) {

  uword pk = XTX.n_rows;
  if (pk == 0) return 1.0;
  if (pk == 1) return as_scalar(XTX(0,0));

  // Lanczos iterations (no reorthogonalization) for the largest eigenvalue of XTX.
  // Unlike the power iteration, it converges fast even when the top eigenvalues are close.
  // Deterministic start: no draw from R's random number generator.
  uword m = std::min(max_it, pk);
  vec v = ones<vec>(pk) + 0.1 * arma::sin(regspace<vec>(1, pk));
  v /= norm(v, 2);
  vec v_old(pk, fill::zeros), alpha(m, fill::zeros), beta(m, fill::zeros);
  double b = 0.0, theta = 0.0, residual = datum::inf;

  for (uword j = 0; j < m; ++j) {
    vec w = XTX * v;
    alpha(j) = dot(w, v);
    w -= alpha(j) * v + b * v_old;
    b = norm(w, 2);
    beta(j) = b;

    bool breakdown = (b <= 1e-12 * std::abs(alpha(j)));
    if (breakdown || j == m - 1 || (j + 1) % 5 == 0) {
      // Largest Ritz value and its residual bound |beta_j * s_j|
      mat T(j + 1, j + 1, fill::zeros);
      T.diag() = alpha.head(j + 1);
      if (j > 0) { T.diag(1) = beta.head(j); T.diag(-1) = beta.head(j); }
      vec ritz; mat S;
      eig_sym(ritz, S, T);
      theta = ritz(j);
      residual = std::abs(b * S(j, j));
      if (breakdown || residual <= tol * theta) break;
    }
    v_old = v;
    v = w / b;
  }

  // theta underestimates the largest eigenvalue; theta + residual bounds it in practice
  return std::max({theta + residual, 1.01 * theta, XTX.diag().max()});
}

uword Optimizer::pgd(
    vec& beta,
    const double& lambda,
    const vec& XTy,
    const mat& XTX,
    std::function<vec(const vec&, double)> proximal_operator,
    const double& accuracy,
    const uword& max_iter,
    const uword m,
    double L_cache) {

  uword p = beta.n_elem;
  mat mat_F(p, m, fill::zeros);
  mat mat_X(p, m, fill::zeros);

  double invL = 1.0 / ((L_cache > 0) ? L_cache : estimate_lipschitz(XTX)); 
  uword iter = 0;
  uword hist = 0; // number of (x_k, f_k) pairs stored since the last restart
  double delta = 2.0 * accuracy;
  double delta_prev = datum::inf;
  bool accelerated = false; // was the current beta obtained by extrapolation?
  vec beta_plain;           // plain proximal gradient step computed at the previous iterate

  while (delta > accuracy && iter < max_iter) {
    // 1. Point fixe standard (G(x))
    vec beta_next = proximal_operator(beta - (XTX * beta - XTy) * invL, lambda * invL);
    vec f_k = beta_next - beta;

    delta = norm(f_k, 2) / invL; // norm of the gradient mapping, invariant to the step size

    // Safeguard: plain proximal gradient steps do not increase the fixed-point residual.
    // If the extrapolated point did, reject it, fall back to the plain step and drop the history.
    if (accelerated && delta > delta_prev) {
      beta = beta_plain ;
      accelerated = false ;
      hist = 0 ;
      delta = delta_prev ;
      iter++;
      continue ;
    }
    delta_prev = delta;
    accelerated = false ;

    if (iter == 0 || m == 0) {
      beta = beta_next;
    } else {
      // 2. Préparation des données pour l'accélération
      uword col_idx = hist % m;       // On stocke l'itéré PRÉCÉDENT
      mat_X.col(col_idx) = beta;      // l'itéré x_k
      mat_F.col(col_idx) = f_k;       // son résidu f_k
      hist++;

      uword current_m = std::min(hist, m);

      if (current_m > 1) {
        // Anderson mixing (type II): dF(:,j) = f_j - f_k (differences from current residual)
        // mat_X.col(j) + mat_F.col(j) = beta_j + (beta_{j+1} - beta_j) = beta_{j+1}
        mat dF = mat_F.cols(0, current_m - 1);
        dF.each_col() -= f_k;
        
        vec gamma;
        if (solve(gamma, dF, -f_k, solve_opts::fast)) {
          vec beta_accel = beta_next;
          for (uword j = 0; j < current_m; ++j) {
            beta_accel += gamma(j) * (mat_X.col(j) + mat_F.col(j) - beta_next);
          }
          beta_plain = std::move(beta_next);
          beta = beta_accel;
          accelerated = true;
        } else {
          beta = beta_next;
        }
      } else {
        beta = beta_next;
      }
    }
    iter++;
    if (iter % 100 == 0) R_CheckUserInterrupt();
  }
  return iter;
}

uword Optimizer::fista(
  vec& beta,
  const double& lambda,
  const vec& XTy,
  const mat& XTX,
  std::function<vec(const vec&, double)> proximal_operator,
  const double& accuracy,
  const uword& max_iter,
  double L_cache) {

  double L = (L_cache > 0) ? L_cache : estimate_lipschitz(XTX);

  vec betak;
  vec betal = beta;
  double delta = 2.0 * accuracy;

  double t0 = 1.0, tk;
  uword iter = 0;
  double invL = 1.0 / L;
  
  while ((delta > accuracy) && (iter < max_iter)) {
    
    // Proximal step
    betak = proximal_operator(betal - (XTX * betal - XTy) * invL, lambda * invL);
    
    // Assess convergence (scaled by L to be invariant to the step size)
    delta = L * norm(beta - betak, 2);
    
    if (dot(betal - betak, betak - beta) > 0) {
      // Adaptive restart (O'Donoghue & Candes, 2015): the momentum points against the
      // gradient mapping, so reset it
      t0 = 1.0;
      betal = betak;
    } else {
      // FISTA update
      tk = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t0 * t0));
      double weight = (t0 - 1.0) / tk;
      
      // Accelerating step
      betal = betak + weight * (betak - beta);
      t0 = tk;
    }
    
    beta = betak;
    iter++;
    
    if (iter % 100 == 0) R_CheckUserInterrupt();
  }
  
  return iter;
}

uword Optimizer::conjugate_gradient(
  vec& x,
  const mat& A,
  const vec& b,
  const double& accuracy,
  const uword& max_iter) {
  
  vec r = b - A * x;
  vec p = r;
  double rs_old = dot(r, r);
  
  if (sqrt(rs_old) < accuracy) return 0;
  
  uword i = 0;
  for (i = 0; i < max_iter; ++i) {
    vec Ap = A * p;
    
    double pAp = dot(p, Ap);
    
    // Handle cases when A is not positive definite
    if (std::abs(pAp) < 1e-16) break;
    
    double alpha = rs_old / pAp;
    
    x += alpha * p;
    r -= alpha * Ap;
    
    double rs_new = dot(r, r);
    
    // Stopping criterion on the residual norm
    if (std::sqrt(rs_new) < accuracy) {
      i++;
      break;
    }
    
    // Update search direction (Fletcher-Reeves)
    p = r + (rs_new / rs_old) * p;
    rs_old = rs_new;
    
    if (i % 100 == 0) R_CheckUserInterrupt();
  }
  
  return i;
}

void Optimizer::optimality_violation(
  const vec& beta,
  const vec& grad,
  const double& lambda,
  const double& gamma,
  const vec& XTy,
  const mat& XTX,
  const double& norm_y,
  uvec A,
  uword type) {
  
  // nu equals the max |gradient|
  double nu   = arma::norm(grad, "inf");
  double loss  = .5 * pow(norm_y, 2) + dot(beta, .5 * XTX * beta - XTy);
  double old_J = J_, old_D = D_;
  J_ = loss - dot(beta, grad(A));
  uvec Ac;
  uword p = grad.n_elem;

  switch (type) {
    case 1: // Grandvalet's bound
      Ac = find(grad > nu);
      D_ = J_ * (1 - lambda/nu) -
        (pow(lambda, 2) / (2*gamma)) * ((lambda*(p - Ac.n_elem))/nu +
        pow(arma::norm(grad(Ac), 2)/nu, 2) - p);
      break;
    case 2: // Fenchel's bound
      if (nu < lambda) nu = lambda;
      D_ = loss * (1 + pow(lambda/nu, 2)) + sum(abs(lambda*beta)) +
        (lambda/nu) * (dot(beta, XTy) - pow(norm_y, 2));
      break;
    default:
      D_ = datum::inf;
      break;
  }

  // keep the smallest bound reached so far for a given lambda value
  if ((old_J < J_) && (old_D - D_) < (old_J - J_)) { D_ = old_D; }
}
