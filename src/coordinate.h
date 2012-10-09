/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _quadrupen_COORDINATE_H
#define _quadrupen_COORDINATE_H

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

List pathwise_enet(vec x0, mat xtx, vec xty, vec xtxw, vec pen, double gam, double eps) ;

#endif
