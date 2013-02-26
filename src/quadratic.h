/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#ifndef _quadrupen_QUADRATIC_H
#define _quadrupen_QUADRATIC_H

#ifndef _RCPP_ARMA_H
#define _RCPP_ARMA_H
#include <RcppArmadillo.h>
#endif

using namespace Rcpp;
using namespace arma;

#include "utils.h"

uvec setdiff(uvec x, uvec y) ;
int quadra_enet(vec& x0, mat& R, vec xty, vec sgn_grd, double &pen, uvec& null) ;
int quadra_breg(vec& beta, const mat& xtx, const vec& xty, double &pen, vec& grd, uvec& B, const int maxit=50) ;

#endif
