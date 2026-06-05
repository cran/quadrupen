/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "PenaltyDense.h"

using namespace Rcpp;
using namespace arma;

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION
template<>
vec DensePenalty<DenseNorm::LINF>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) % w);
}

template<>
vec DensePenalty<DenseNorm::LINF>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::abs(x) / w);
}

template<>
double DensePenalty<DenseNorm::LINF>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(max(elt_norm(x, w)));
}

template<>
double DensePenalty<DenseNorm::LINF>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_dual_norm(x,w))) ;
}

template<>
vec DensePenalty<DenseNorm::LINF>::proximal(const vec& x, double lambda, const vec& w) {
  uword p = x.n_elem;
  vec abs_x = arma::abs(x);

  // If x is already in the dual ball {y: sum(|y_i|/w_i) <= lambda}, proximal is 0
  if (accu(abs_x / w) <= lambda) {
    return zeros<vec>(p);
  }

  // Project x onto {y : sum(|y_i|/w_i) <= lambda} via KKT conditions:
  //   y_i = sign(x_i) * (|x_i| - mu/w_i)_+
  // where mu satisfies sum((|x_i|/w_i - mu/w_i^2)_+) = lambda.
  // Breakpoints are at mu = |x_i|*w_i; sort these descending.
  uvec ord = sort_index(abs_x % w, "descend");

  double cumsum_z   = 0.0;  // sum of |x_i|/w_i for the active set
  double cumsum_iw2 = 0.0;  // sum of 1/w_i^2 for the active set
  double mu = 0.0;

  for (uword j = 0; j < p; ++j) {
    uword k = ord(j);
    cumsum_z   += abs_x(k) / w(k);
    cumsum_iw2 += 1.0 / (w(k) * w(k));
    double t = (cumsum_z - lambda) / cumsum_iw2;

    if (j < p - 1 && abs_x(ord(j + 1)) * w(ord(j + 1)) <= t) {
      mu = t;
      break;
    }
    if (j == p - 1) mu = t;
  }

  // prox_i = x_i - y_i = sign(x_i) * min(|x_i|, mu/w_i)
  return sign(x) % arma::min(abs_x, mu / w);
}

// ______________________________________________________
// L2 NORM SQUARED A.K.A RIDGE
template<>
vec DensePenalty<DenseNorm::L2>::elt_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) % w);
}

template<>
vec DensePenalty<DenseNorm::L2>::elt_dual_norm(const vec& x, const vec& w, double lambda) {
  return(arma::pow(x, 2) / w);
}

template<>
double DensePenalty<DenseNorm::L2>::pen_norm(const vec& x, const vec& w, double lambda) {
  return(sum(elt_norm(x, w)));
}

// Slight abuse: take the Fenchel conjugate of L2^2
template<>
double DensePenalty<DenseNorm::L2>::dual_norm(const vec& x, const vec& w, double lambda) {
  return(.25*sum(elt_norm(x, w))) ;
}

template<>
vec DensePenalty<DenseNorm::L2>::proximal(const vec& x, double lambda, const vec& w) {
  return(x / (1+2*lambda*w));
}
