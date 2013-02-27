/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "quadratic.h"

int quadra_enet(vec &x0,
		mat &R,
		mat &xAtxA,
		vec  xty,
		vec  sgn_grd,
		double &pen ,
		uvec   &null,
		bool   usechol,
		double tol) {

  uword iter = 1; // current iterate
  
  uvec A = find(abs(x0) > ZERO) ; // vector of active variables
  vec theta = -sgn_grd   ; // vector of sign of the solution
  theta.elem(A)   = sign(x0.elem(A));
  
  // Solving the quadratic problem
  vec x1 ;
  if (usechol) {
    x1 = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen * theta));
  } else {
    x1 = cg(xAtxA, xty - pen*theta, x0, tol) ;
  }
  
  // Check for swapping variables
  uvec swap = find(abs(sign(x1.elem(A)) - theta.elem(A)) > ZERO);
  if (swap.is_empty()) {
    null = swap ; // this is empty
    x0 = x1;
  } else {
    swap = A.elem(swap);
    colvec x0_swap = x0.elem(swap);
    colvec x1_swap = x1.elem(swap);
    // first, go to zero for the swapped variable which cost the minimum
    vec gamma = -x0_swap / (x1_swap-x0_swap);
    uword i_swap ;
    double scale = gamma.min(i_swap) ;
    null = swap[i_swap];
    x1 = x0 + (x1-x0) * scale ;
    // second, solve the problem after swaping the signs of the
    // increminated variable
    x0 = x1;
    x0(null[0]) = -x1_swap[i_swap];

    A = find(abs(x0) > ZERO) ; // vector of active variables
    theta = -sgn_grd        ; // vector of sign of the solution
    theta.elem(A)   = sign(x0.elem(A));

    vec x2 ;
    if (usechol) {
      x2 = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen * theta));
    } else {
      x2 = cg(xAtxA, xty - pen*theta, x1, tol) ;
    }
    iter++;

    // This is the gradient on the active part of the parameters
    vec grd = -xty + xAtxA * x2;
    // if the sign is coherent, keep that one...
    if (fabs(grd(null[0]) + pen * as_scalar(sign(x2(null)))) <= ZERO) {
      null = swap; // this is empty
      x0 = x2 ;
    } else {
      // otherwise, backtrack to x1
      x0 = x1 ;
    }
  }
  return(iter);
}

int quadra_breg(vec    &beta,
		const mat& xtx,
		const vec& xty,
		double &pen   ,
		vec    &grd   ,
		uvec   &B     ,
		const int maxit) {

  const double zero = 2e-16     ;
  int p        = beta.n_elem    ; // size of the problem
  int iter     = 0              ; // count the number of systems solved
  double bound ; //
  uvec all(p)           ;
  for (int i=0;i<p;i++){all(i) = i;}
  uvec I             ; // guys living in between the supremum
  uvec toB           ; // guys reaching the boundary after optimization
  uvec toI           ; // guys leaving the boundary after optimization
  vec  theta = -sign(grd.elem(B)) ;

  vec XX_B   ;
  mat XX     ;
  mat XX_II  ;
  mat R      ;
  vec b      ;
  vec tmp    ;

  I = setdiff(all,B);
  toB = B;

  while (toB.n_elem > 0 & iter < maxit) {

    iter++;
    //
    // SOLVE THE QUADRATIC PROBLEM
    //

    // Constructing the system (KKT)
    XX_B = xtx.cols(B) * theta;

    if (I.is_empty()) {
      XX = sum(theta % XX_B.elem(B),0);
      b  = (dot(theta, xty.elem(B))-pen);
      beta.elem(B) = theta * (b/XX) ;
      // keep a trace of the current boundary
    } else {
      if (I.n_elem > 1) {
        XX_II = xtx.submat(I,I) ;
      } else {
        XX_II = xtx(I,I) ;
      }
      XX   = join_rows(
             join_cols(sum(theta % XX_B.elem(B),0),XX_B.elem(I)),
             join_cols(strans(XX_B.elem(I)), XX_II));

      vec b = zeros<vec>(I.n_elem + 1) ;
      b[0] = dot(theta, xty.elem(B))-pen;
      b.subvec(1,b.n_elem-1) = xty.elem(I) ;

      // Solving with Cholesky factorization...
      R = chol(XX) ;
      tmp = solve(trimatu(R), solve(trimatl(strans(R)),b)) ;
      beta.elem(B) = theta * tmp[0] ;
      beta.elem(I) = tmp.subvec(1,tmp.n_elem-1) ;
    }
    // keep a trace of the current boundary
    bound = max(abs(beta.elem(B)));
    //
    // VARIABLES REACHING THE BOUNDARY
    //
    toB = find(abs(beta) > bound);
    beta.elem(toB) = bound * sign(beta.elem(toB));
    B = unique(join_cols(B,toB));
    I = setdiff(all,B);
    theta = sign(beta.elem(B)); // sign of the guys reaching the supremum
  }

  grd = -xty + xtx * beta ;

  //
  // VARIABLE LEAVING THE BOUDARY
  //
  toI = find(abs(theta + sign(grd.elem(B))) > zero);
  if (!toI.is_empty()) {
    toI = B.elem(toI);
  }
  if (!toI.is_empty()) {
    B = setdiff(B,toI);
     I = setdiff(all,B);
  }
  // If everyone is leaving the boundary, that's a deal...
  if (B.is_empty()) {
    throw std::runtime_error("Too much unstability");
  }

  return(iter) ;
}

uvec setdiff(uvec x, uvec y) {

  uword ind_x = 0;
  uword ind_y = 0;
  uword ind_z = 0;
  uword end_x = x.n_elem;
  uword end_y = y.n_elem;
  uvec z = zeros<uvec>(x.n_elem-y.n_elem);

  if (y.is_empty()) {
    z = x;
  } else {
    while (ind_y != end_y) {
      if ( x(ind_x) < y(ind_y) ) {
	z(ind_z) = x(ind_x);
	ind_z++;
	ind_x++;
      } else if ( y(ind_y) < x(ind_x) ) {
	ind_y++;
      } else {
	ind_x++;
	ind_y++;
      }
      R_CheckUserInterrupt();
    }
    while(ind_x != end_x) {
      z(ind_z) = x(ind_x);
      ind_z++;
      ind_x++;
    }
  }
  return(z);
}
