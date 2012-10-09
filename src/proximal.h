/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_PROXIMAL_H
#define _quadrupen_PROXIMAL_H

#include <RcppArmadillo.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

List fista_lasso(vec x0, mat xtx, vec xty, vec pen, double L0, double eps) ;

#endif
