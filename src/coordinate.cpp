/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "coordinate.h"

List pathwise_enet(vec  x0,
		   mat xtx,
		   vec xty,
		   vec xtxw,
		   vec pen,
		   double gam   ,
		   double eps    ) {
  
  colvec xk = x0 ; // output vector
  int j, i     = 0     ; // current iterate
  int max_iter = 10000 ; // max. number of iteration
  double delta = 2*eps ; // change in beta
  double u, d          ; // temporary scalar

  while ((delta > eps/x0.n_elem ) && (i < max_iter)) {

    delta = 0;
    for (j=0; j<x0.n_elem; j++) {
      // Soft thresholding operator
      u = x0(j) * (1+gam) + xty(j) - xtxw(j) ;
      xk(j)  = fmax(1-pen(j)/fabs(u),0) * u/(1+gam) ;
      d = xk(j)-x0(j);
      delta += pow(d,2);
      xtxw  += d*xtx.col(j) ;
    }
    
    // preparing next iterate
    delta = sqrt(delta);
    x0 = xk;
    i++;

    R_CheckUserInterrupt();    
  }

 return List::create(Named("iter")  = i,
			    Named("xk")   = xk,
			    Named("xtxw") = xtxw) ;
}
