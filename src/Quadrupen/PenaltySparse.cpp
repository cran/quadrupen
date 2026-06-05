/*
* Author: Julien CHIQUET
*         MIA Paris-Saclay
*/
  
#include "PenaltySparse.h"
#include "PenaltyUtils.h"

using namespace Rcpp;
using namespace arma;

// ─────────────────────────────────────────────
// Shared Helpers 
// ─────────────────────────────────────────────
namespace {

// Violation KKT d'une variable inactive : identique L1/SCAD/MCP au voisinage de 0
inline double violation_inactive(double gj, double wj, double lambda) {
  return std::abs(gj) / wj - lambda;  // > 0  ↔  violation
}

// Violation KKT d'une variable active pour une pénalité quelconque :
// |g_j + d_j * sign(β_j)|  où d_j est la dérivée effective de la pénalité
inline double violation_active(double gj, double bj, double dj) {
  double sj = (bj >= 0.0) ? 1.0 : -1.0;
  return std::abs(gj + dj * sj);
}

} // namespace

// ______________________________________________________
// L1 NORM A.K.A LASSO
template<>
vec SparsePenalty<SparseNorm::L1>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec SparsePenalty<SparseNorm::L1>::proximal(const vec& x, double lambda, const vec& w) {
  return soft_threshold(x, lambda * w) ;
}

template<>
vec SparsePenalty<SparseNorm::L1>::derivative(const vec& beta, double lambda, const vec& w) {
  return lambda * w;
}

template<>
vec SparsePenalty<SparseNorm::L1>::optimality(
    const vec& grad, double lambda, const vec& w,
    const vec& beta, const uvec& A) const {
  
  // Condition par défaut : violation L1 pour toutes les variables
  vec viol = arma::abs(grad) / w - lambda;
  
  // Correction pour les variables actives
  // KKT exacte : g_j + λ w_j sign(β_j) = 0
  if (!beta.is_empty() && !A.is_empty()) {
    for (uword k = 0; k < A.n_elem; ++k) {
      double bj = beta[k];
      if (std::abs(bj) > 0.0) {
        double dj = lambda * w[A[k]];  // dérivée L1
        viol[A[k]] = violation_active(grad[A[k]], bj, dj);
      }
    }
  }
  return viol;
}

// ______________________________________________________
// MCP

template<>
vec SparsePenalty<SparseNorm::MCP>::elt_norm(const vec& x, const vec& w, double lambda) {
  vec res = zeros<vec>(x.n_elem);
  for(uword i=0; i<x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    if (abs_xi <= eta_ * lambda) {
      res[i] = lambda * w[i] * abs_xi - (abs_xi * abs_xi) / (2.0 * eta_);
    } else {
      res[i] = 0.5 * eta_ * lambda * lambda * w[i]; // constant value (plateau)
    }
  }
  return res;
}

template<>
vec SparsePenalty<SparseNorm::MCP>::proximal(const vec& x, double lambda, const vec& w) {
  vec res = zeros<vec>(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      res[i] = 0.0;
    } else if (abs_xi <= eta_ * l) {
      // equivalent to / (1 - 1/eta_), but avoids singularity at eta_=1
      res[i] = std::copysign((abs_xi - l) * eta_ / (eta_ - 1.0), x[i]);
    } else {
      res[i] = x[i];
    }
  }
  return res;
}

template<>
vec SparsePenalty<SparseNorm::MCP>::derivative(const vec& beta, double lambda, const vec& w) {
  vec d = zeros<vec>(beta.n_elem);
  for(uword i=0; i<beta.n_elem; ++i) {
    double abs_b = std::abs(beta[i]);
    if (abs_b < eta_ * lambda * w[i]) {
      d[i] = lambda * w[i] - abs_b / eta_;
    } else {
      d[i] = 0.0;
    }
  }
  return d;
}

// ─────────────────────────────────────────────
// MCP  (paramètre γ = eta_ > 1)
//
// Dérivée effective selon |β_j| :
//   |β| ≤ γλw  →  d = λw − |β|/γ   (zone linéaire décroissante)
//   |β| > γλw  →  d = 0             (plateau)
// ─────────────────────────────────────────────
template<>
vec SparsePenalty<SparseNorm::MCP>::optimality(
    const vec& grad, double lambda, const vec& w,
    const vec& beta, const uvec& A) const {
  
  vec viol = arma::abs(grad) / w - lambda;
  
  if (!beta.is_empty() && !A.is_empty()) {
    for (uword k = 0; k < A.n_elem; ++k) {
      uword j   = A[k];
      double bj = beta[k];
      if (std::abs(bj) < 1e-10) continue;
      
      double abs_bj = std::abs(bj);
      double lj     = lambda * w[j];
      double dj;
      
      if (abs_bj <= eta_ * lj) {
        dj = lj - abs_bj / eta_;   // zone MCP
      } else {
        dj = 0.0;                  // plateau
      }
      
      viol[j] = violation_active(grad[j], bj, dj);
    }
  }
  return viol;
}
// ______________________________________________________
// SCAD

template<>
vec SparsePenalty<SparseNorm::SCAD>::elt_norm(const vec& x, const vec& w, double lambda) {
  vec res = zeros<vec>(x.n_elem);
  for(uword i=0; i<x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      res[i] = l * abs_xi;
    } else if (abs_xi <= eta_ * l) {
      res[i] = (2.0 * eta_ * l * abs_xi - abs_xi * abs_xi - l * l) / (2.0 * (eta_ - 1.0));
    } else {
      res[i] = (l * l * (eta_ + 1.0)) / 2.0;
    }
  }
  return res;
}

template<>
vec SparsePenalty<SparseNorm::SCAD>::proximal(const vec& x, double lambda, const vec& w) {
  vec res = zeros<vec>(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) {
    double abs_xi = std::abs(x[i]);
    double l = lambda * w[i];
    if (abs_xi <= l) {
      // Case 1 : Sparse zone (Lasso)
      res[i] = 0.0;
    } else if (abs_xi <= 2.0 * l) {
      // Case 2 : Soft-thresholding standard
      res[i] = std::copysign(abs_xi - l, x[i]);
    } else if (abs_xi <= eta_ * l) {
      // Case 3 : Transition (SCAD Threshoholding
      // ((eta - 1) * x - sign(x) * eta * l) / (eta - 2)
      double val = ((eta_ - 1.0) * abs_xi - eta_ * l) / (eta_ - 2.0);
      res[i] = std::copysign(std::max(0.0, val), x[i]);
    } else {
      // Case 4 : No bias
      res[i] = x[i];
    }
  }
  return res;  
}

template<>
vec SparsePenalty<SparseNorm::SCAD>::derivative(const vec& beta, double lambda, const vec& w) {
  vec d = zeros<vec>(beta.n_elem);
  for(uword i=0; i<beta.n_elem; ++i) {
    double abs_b = std::abs(beta[i]);
    double l = lambda * w[i];
    if (abs_b <= l) {
      d[i] = l;
    } else if (abs_b <= eta_ * l) {
      d[i] = (eta_ * l - abs_b) / (eta_ - 1.0);
    } else {
      d[i] = 0.0;
    }
  }
  return d;
}

// ─────────────────────────────────────────────
// SCAD  (paramètre a = eta_ > 2, typiquement 3.7)
//
// Dérivée effective selon |β_j| :
//   |β| ≤ λw        →  d = λw           (zone L1)
//   λw < |β| ≤ aλw  →  d = (aλw−|β|)/(a−1)  (zone concave)
//   |β| > aλw       →  d = 0            (plateau : pas de pénalité)
// ─────────────────────────────────────────────
template<>
vec SparsePenalty<SparseNorm::SCAD>::optimality(
    const vec& grad, double lambda, const vec& w,
    const vec& beta, const uvec& A) const {
  
  vec viol = arma::abs(grad) / w - lambda;  // condition inactive (L1 en 0)
  
  if (!beta.is_empty() && !A.is_empty()) {
    for (uword k = 0; k < A.n_elem; ++k) {
      uword j   = A[k];
      double bj = beta[k];
      if (std::abs(bj) < 1e-10) continue; // considered inactive
      
      double abs_bj = std::abs(bj);
      double lj     = lambda * w[j];
      double dj;
      
      if (abs_bj <= lj) {
        dj = lj;                                   // zone L1
      } else if (abs_bj <= eta_ * lj) {
        dj = (eta_ * lj - abs_bj) / (eta_ - 1.0);  // zone concave
      } else {
        dj = 0.0;                                  // plateau
      }
      
      viol[j] = violation_active(grad[j], bj, dj);
    }
  }
  return viol;
}
