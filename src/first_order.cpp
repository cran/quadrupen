/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "first_order.h"

int fista_lasso(vec  &x0   ,
		mat  &xtx  ,
		vec xty  ,
		double &pen,
		uvec   &null,
		double &L0 ,
		double eps) {

  colvec xk = x0  ; // output vector
  colvec  s = x0  ;
  int iter = 0, j = 0      ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration
  double L  = L0        ; // Lipschitz's constant

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;
  colvec df ;

  double l_num, l_den ;

  while ((delta > eps/x0.n_elem ) && (iter < max_iter)) {
    
    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    df = - xty + xtx * s ;
    
    // Line search over L
    while(!found) {
      // apply proximal operator of the Lasso
      xk = s - df/L ;
      for (j=0; j<x0.n_elem; j++) {
	xk(j) = fmax(0,1-(pen/L)/fabs(xk(j))) * xk(j);
      }

      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(df, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));

      if ((L * l_den >= l_num) || (sqrt(l_den) < eps)) {
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
    iter++;

    R_CheckUserInterrupt();
  }

  null = sort(find(abs(xk) + (abs(df) - pen) < ZERO),1) ;

  return(iter);
}

int fista_breg(vec    &x0,  
	       const mat &xtx,   
	       const vec &xty,   
	       vec    &grd,
	       double &pen,
	       double &L0 ,  
	       double eps) {
  
  colvec xk = x0        ; // output vector
  colvec  s = x0        ;
  int iter = 0          ; // current iterate
  double delta = 2*eps  ; // change in beta
  int max_iter  = 10000 ; // max. number of iteration

  double t0 = 1.0, tk;
  bool found=false;
  double f0, fk ;

  double l_num, l_den ;

  while ((delta > eps*eps) && (iter < max_iter)) {

    f0 = as_scalar(.5 * strans(s) * xtx * s - strans(xty) * s) ;
    grd = - xty + xtx * s ;
    
    // Line search over L
    while(!found) {
      xk = proximal_inf(s - grd/L0, pen/L0);
      
      fk = as_scalar(.5 * strans(xk) * xtx * xk - strans(xty) * xk) ;
      l_num = as_scalar(2 * (fk - f0 - dot(grd, xk-s) ));
      l_den = as_scalar(pow(norm(xk-s,2),2));
      
      if ((L0 * l_den >= l_num) || (sqrt(l_den) < eps)) {
	found = true;
      } else {
	L0 = fmax(2*L0, l_num/l_den);
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
    iter++;

    R_CheckUserInterrupt();
  }

  return(iter) ;
}

vec proximal_inf(vec v,
		 double lambda) {
  
  int p = v.n_elem;
  vec u, proj;
  vec out = zeros<vec>(p);
  
  if ( as_scalar(sum(abs(v) / lambda)) >= 1) {   
    
    // Reordonnons les valeurs absolues
    u = sort(abs(v),1);
    
    // valeurs des coordonnées projetees si non nulles (problème dual)
    proj = (cumsum(u) - lambda)/linspace<vec>(1,p,p);
    
    // selection des coordonnees non nulles (problème dual)
    uvec maxs = sort(find(u-proj>ZERO),1) ;
    double thresh = proj[maxs[0]];
    
    // solution du programme primal
    // On garde les valeurs les plus petites, et on seuille le reste à une valeur commune +- thresh
    for (int k=0; k < p;k++) {
      if (fabs(v(k)) > ZERO) {
	if (v(k) > 0) {
	  out(k) = fmin(fabs(v(k)),thresh);
	} else {
	  out(k) = -fmin(fabs(v(k)),thresh);
	}
      }
    }
  }
  return(out);
}

int pathwise_enet(vec&  x0,
		  mat& xtx,
		  vec xty,
		  vec& xtxw,
		  double& pen,
		  uvec &null,
		  double& gam   ,
		  double eps    ) {
  
  colvec xk = x0 ; // output vector
  int j, iter  = 0     ; // current iterate
  int max_iter = 10000 ; // max. number of iteration
  double delta = 2*eps ; // change in beta
  double u, d          ; // temporary scalar

  while ((delta > eps/x0.n_elem ) && (iter < max_iter)) {

    delta = 0;
    for (j=0; j<x0.n_elem; j++) {
      // Soft thresholding operator
      u = x0(j) * (1+gam) + xty(j) - xtxw(j) ;
      xk(j)  = fmax(1-pen/fabs(u),0) * u/(1+gam) ;
      d = xk(j)-x0(j);
      delta += pow(d,2);
      xtxw  += d*xtx.col(j) ;
    }
    
    // preparing next iterate
    delta = sqrt(delta);
    x0 = xk;
    iter++;

    R_CheckUserInterrupt();    
  }

  null = sort(find(abs(xk) + (abs(-xty + xtxw) - pen) < ZERO),1) ;
  return(iter);
}
