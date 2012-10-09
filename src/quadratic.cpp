/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#include "quadratic.h"

List quadra_enet(vec x0,
		      mat  R,
		      vec  xty,
		      vec  sgn_grd,
		      vec  pen,
		      double eps    ) {

  uword iter = 0; // current iterate

  uvec A = find(abs(x0) > eps) ; // vector of active variables
  vec theta = -sgn_grd   ; // vector of sign of the solution
  theta.elem(A)   = signs(x0.elem(A),eps);
  
  // Solving the quadratic problem
  vec x1  = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen % theta));
  
  // Check for swapping variables
  uvec swap = find(abs(signs(x1.elem(A), eps) - theta.elem(A)) > eps);
  if (swap.is_empty()) {
    return List::create(iter, x1) ;
  } else {
    swap = A.elem(swap);
    colvec x0_swap = x0.elem(swap);
    colvec x1_swap = x1.elem(swap);
    // first, go to zero for the swapped variable which cost the minimum
    vec gamma = -x0_swap / (x1_swap-x0_swap);
    uword i_swap ;
    double scale = gamma.min(i_swap) ;
    uword swap_min = swap[i_swap];
    x1 = x0 + (x1-x0) * scale ;
    // second, solve the problem after swaping the signs of the
    // increminated variable
    x0 = x1;
    x0(swap_min) = -x1_swap[i_swap];
    
    A = find(abs(x0) > eps) ; // vector of active variables
    theta = -sgn_grd        ; // vector of sign of the solution
    theta.elem(A)   = signs(x0.elem(A),eps);

    iter++;
    vec x2 = solve(trimatu(R), solve( trimatl(strans(R)), xty - pen % theta));
    return List::create(iter, x1, x2, swap_min);    
  }  
}
