/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#ifndef _quadrupen_QUADRATIC_H
#define _quadrupen_QUADRATIC_H

#include <RcppArmadillo.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

uvec setdiff(uvec x, uvec y) ;
int quadra_enet(vec& x0, mat& R,  mat& xAtxA, vec xty, vec sgn_grd, double &pen, uvec& null, bool usechol, double tol) ;
int quadra_breg(vec& beta, const mat& xtx, const vec& xty, double &pen, vec& grd, uvec& B, const int maxit=50) ;

#endif
