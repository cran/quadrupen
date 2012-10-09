/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "proximal.h"

List fista_lasso(vec  x0,
		 mat xtx,
		 vec xty,
		 vec pen,
		 double L0    ,
		 double eps   ) {

  colvec xk = x0  ; // output vector
  colvec  s = x0  ;
  int i = 0, j = 0      ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration
  double L  = L0        ; // Lipschitz's constant
  // double L  = eig_sym(xtx).max() ; // Lipschitz's constant

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  colvec df ;

  double l_num, l_den ;

  while ((delta > eps/x0.n_elem ) && (i < max_iter)) {

    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    df = (xtx * s - xty) ;
    
    // Line search over L
    while(!found) {
      // apply proximal operator of the Lasso
      xk = s - df/L ;
      for (j=0; j<x0.n_elem; j++) {
	xk(j) = fmax(0,1-(pen(j)/L)/fabs(xk(j))) * xk(j);
      }

      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(df, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));

      if (L * l_den >= l_num  | sqrt(l_den) < eps) {
	found = true;
      } else {
	L = fmax(2*L, l_num/l_den);
      }
      
      R_CheckUserInterrupt();
    }
    
    // updating t 
    tk = 0.5 * (1+sqrt(1+4*t0*t0));
    
    // updating s
    s = xk + (t0-1)/tk * ( xk - x0 );
    
    // preparing next iterate
    delta = sqrt(l_num);
    x0 = xk;
    t0 = tk;
    found = false;
    i++;

    R_CheckUserInterrupt();
  }

  return List::create(i, xk, L) ;
}
