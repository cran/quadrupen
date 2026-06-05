/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

// Shared proximal utility functions used across penalty implementations.

#pragma once
#include <armadillo>

// Signed soft thresholding: sign(x) * max(|x| - t, 0), element-wise.

// Scalar threshold (uniform across all elements).
inline arma::vec soft_threshold(const arma::vec& x, double t) {
  return arma::sign(x) % arma::clamp(arma::abs(x) - t, 0.0, arma::datum::inf);
}

// Per-element threshold.
inline arma::vec soft_threshold(const arma::vec& x, const arma::vec& t) {
  return arma::sign(x) % arma::clamp(arma::abs(x) - t, 0.0, arma::datum::inf);
}
