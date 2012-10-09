/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_UTILS_H
#define _quadrupen_UTILS_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

vec signs(vec  x, double zero) ;

double sign(double  x, double zero) ;

mat cholupdate(mat R , mat XtX) ;

mat choldowndate(mat R, int j) ;

#endif 
