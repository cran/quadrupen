/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_UTILS_H
#define _quadrupen_UTILS_H

#ifndef _RCPP_ARMA_H
#define _RCPP_ARMA_H
#include <RcppArmadillo.h>
#endif 

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

inline double sign(double x) {return((x > ZERO) ? 1.0 : ((-x > ZERO) ? -1.0 : 0.0)) ;}
inline vec signs(vec x) {
  vec signs = zeros<vec>(x.n_elem);
  for (int j=0; j<x.n_elem; j++) {
    signs(j) = (x(j) > ZERO) ? 1.0 : ((-x(j) > ZERO) ? -1.0 : 0.0);
  }
  return(signs);
}
#endif 
