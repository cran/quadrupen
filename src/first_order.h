/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_FIRSTORDER_H
#define _quadrupen_FIRSTORDER_H

#ifndef _RCPP_ARMA_H
#define _RCPP_ARMA_H
#include <RcppArmadillo.h>
#endif 

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

int pathwise_enet(vec& x0, mat& xtx, vec xty, vec& xtxw, double& pen, uvec &null, double& gam, double eps) ;
int fista_lasso(vec &x0, mat &xtx, vec xty, double &pen, uvec &null, double &L0, double eps) ;
int fista_breg (vec &x0, const mat &xtx, const vec &xty, vec& grd, double &pen, double &L0, double eps) ;
vec proximal_inf(vec v, double lambda) ;

#endif
