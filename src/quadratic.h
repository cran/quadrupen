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

List quadra_enet(vec x0, mat R, vec xty, vec sgn_grd, vec pen, double eps) ; 

#endif
