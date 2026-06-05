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
  
}

double Optimizer::estimate_lipschitz(
  const mat& XTX,
  uword max_it,
  double tol) {

  uword pk = XTX.n_rows;
  if (pk == 0) return 1.0;
  if (pk == 1) return as_scalar(XTX(0,0));

  // Warm-start from previous eigenvector when size matches; random init otherwise
  vec q;
  if (q_lipschitz_.n_elem == pk) {
    q = q_lipschitz_;
  } else {
    q = randu<vec>(pk);
    q /= norm(q, 2);
  }

  double lambda = 0.0;
  double lambda_old = 0.0;

  for (uword i = 0; i < max_it; ++i) {
    vec z = XTX * q;

    // Largest eigenvalue (simplified Rayleigh quotient since ||q||=1)
    lambda = dot(q, z);
    if (i > 0 && std::abs(lambda - lambda_old) < tol * lambda) {
      break;
    }
    lambda_old = lambda;

    double n = norm(z, 2);
    if (n > 1e-15) {
      q = z / n;
    } else {
      break;
    }
  }

  q_lipschitz_ = q; // save for next call

  // Safety margin for 1/L
  return lambda * 1.01;
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
  double delta = 2.0 * accuracy;
  
  while (delta > accuracy && iter < max_iter) {
    // 1. Point fixe standard (G(x))
    vec beta_next = proximal_operator(beta - (XTX * beta - XTy) * invL, lambda * invL);
    vec f_k = beta_next - beta;
    
    delta = norm(f_k, 2);
    
    if (iter == 0 || m == 0) {
      beta = beta_next;
    } else {
      // 2. Préparation des données pour l'accélération
      uword col_idx = (iter - 1) % m; // On stocke l'itéré PRÉCÉDENT
      mat_X.col(col_idx) = beta;      // l'itéré x_k
      mat_F.col(col_idx) = f_k;       // son résidu f_k
      
      uword current_m = std::min(iter, m);

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
          beta = beta_accel;
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
    
    // FISTA update
    tk = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t0 * t0));
    double weight = (t0 - 1.0) / tk;
    
    // Accelerating step
    betal = betak + weight * (betak - beta);
    
    // Assess convergence
    delta = norm(beta - betak, 2);
    
    beta = betak;
    t0 = tk;
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
