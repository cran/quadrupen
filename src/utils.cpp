/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "utils.h"

vec signs(vec  x,
		double zero) {
  
  vec signs = zeros<vec>(x.n_elem);
  
  for (int j=0; j<x.n_elem; j++) {
    if (x(j) > zero) {
      signs(j) = 1;
    } else if (x(j) < -zero) {
      signs(j) = -1;      
    }
  }
  return(signs);
}

double sign(double  x,
	    double zero) {
  
  double sign = 0.0;

  if (x > zero) {
    sign = 1;
  } else if (x < -zero) {
    sign = -1;      
  }
  return(sign);
}

mat cholupdate(mat R , mat XtX) {
  int p = XtX.n_cols;

  if (p == 1) {
    R = sqrt(XtX);
  } else {
    colvec rp  = zeros<colvec>(p,1);
    rp.subvec(0,p-2) = solve (trimatl(strans(R)), XtX.submat(0,p-1,p-2,p-1));
    rp(p-1) = sqrt(XtX(p-1,p-1) - dot(rp,rp));
    R = join_rows( join_cols(R, zeros<mat>(1,p-1)) , rp);    
  }

  return(R);
  
}

mat choldowndate(mat R, int j) {

  vec x = zeros<vec>(2,1);
  mat G = zeros<mat>(2,2);
  
  R.shed_col(j);
  int p = R.n_cols;
  double r;
  for (int k=j; k<p; k++) {
    x = R.submat(k,k,k+1,k);
    
    if (x[1] != 0) {
      r = norm(x,2);
      G <<  x(0) << x(1) << endr
	<< -x(1) << x(0) << endr;
      G = G / r;
      x(0) = r; x(1) = 0;
    } else {
      G = eye(2,2);
    }
    R.submat(k,k,k+1,k) = x;
    if (k < p-1) {
      R.submat(k,k+1,k+1,p-1) = G * R.submat(k,k+1,k+1,p-1);
    }
  }
  R.shed_row(p);
  
  return (R);
  
}
