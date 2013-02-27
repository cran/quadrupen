/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "bounded_reg.h"

using namespace Rcpp;
using namespace arma;

SEXP bounded_reg(SEXP X, SEXP XTY, SEXP STRUCT, SEXP LAMBDA1, SEXP LAMBDA2, SEXP XBAR, SEXP NORMX, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP BULLETPROOF) {

  // disable messages being printed to the err2 stream (armadillo's runtime error)
  std::ostream nullstream(0);
  set_stream_err2(nullstream);

  // Reading input variables
  mat    x        = as<mat>    (X)           ; // densely encoded SCALED design matrix
  vec    lambda1  = as<vec>    (LAMBDA1)     ; // penalty levels
  double lambda2  = as<double> (LAMBDA2)     ; // penalty levels
  mat    struc    = as<mat>    (STRUCT)      ; // the structuring matrix
  vec    xbar     = as<vec>    (XBAR)        ; // mean of the predictors
  vec    normx    = as<vec>    (NORMX)       ; // norm of the predictors
  vec    weights  = as<vec>    (WEIGHTS)     ; // norm of the predictors
  bool   naive    = as<bool>   (NAIVE)       ; // naive or not (renormalize coefficients?)
  vec    xty      = as<vec>    (XTY)         ; // reponses to predictors vector
  double eps      = as<double> (EPS)         ; // precision required
  uword  fun      = as<int>    (FUN)         ; // solver (0=quadra, 1=fista)
  int    verbose  = as<int>    (VERBOSE)     ; // int for verbose mode (0/1/2)
  bool   bullet   = as<bool>   (BULLETPROOF) ; // int for verbose mode (0/1/2)
  uword  max_iter = as<int>    (MAXITER)     ; // max # of iterates of the active set
  uword  max_feat = as<int>    (MAXFEAT)     ; // max # of variables activated

  // Managing dense coding for the design matrix (sparsity is useless here)
  x = x - sqrt(weights) * trans(xbar)        ; // dense encoding of the design matrix
  mat xtx = strans(x) * x + lambda2*struc    ; // t(x) * x + lambda2.S

  // Initializing "first level" variables (outside of the lambda1 loop)
  uword  p        = xty.n_elem             ; // problem size
  uword  n_lambda = lambda1.n_elem         ; // # of penalty
  colvec beta     = zeros<vec>(p)          ; // vector of current parameters
  uvec   B          (p)                    ; // guys reaching the boundary
  for (int i=0;i<p;i++){B(i) = i;}
  mat    coef                              ; // matrice of solution path
  vec    grd                               ; // smooth part of the gradient
  vec    max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
  vec    converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
  uvec   it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda1
  uvec   it_optim                          ; // # of loop in the optimization process for each loop of the active se
  bool  success_optim = true               ; // was the internal system resolution successful?
  vec    timing      (n_lambda)            ; // succesive timing in
  wall_clock timer                         ; // clock
  mat   iB                                 ; // contains row indices of the bounded variables
  mat   jB                                 ; // contains column indices of the bounded variables

  // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
  int  nbr_opt  = 0                        ; // # of current calls to the optimization routine
  double L      = max (xtx.diag())         ; // Lipschitz coefficients

  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (int m=0; m<n_lambda; m++) {
    if (verbose == 2) {Rprintf("\n lambda1 = %f",lambda1(m)) ;}

    // _____________________________________________________________
    //
    // START THE ALGORITHM
    // _____________________________________________________________
    //

    // smooth part of the gradient
    grd     = -xty + xtx * beta        ;
    // dual norm of the gradient
    max_grd[m] = as_scalar(sum(abs(grd)) - lambda1[m]);
    if (max_grd[m] < 0) {
      max_grd[m] = 0;
    }
    
    while ((max_grd[m] > eps) && (it_active[m] <= max_iter)) {
      // _____________________________________________________________
      //
      // (1) KKT/SYSTEM RESOLUTION
      // _____________________________________________________________

      // save the number of iterates performed along the optimization process
      it_optim.reshape(nbr_opt + 1,1) ;

      switch (fun) {
      case 1 :
	it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda1[m], L, pow(eps,2));
        // Evaluating the set of variable reaching the boundary
	B    = find(abs(abs(beta) - max(abs(beta))) < ZERO );
	break;
      default:
	try {
	  // If no convergence up to now...
	  if (it_active[m] == max_iter) {
	    throw std::runtime_error("Fail to converge...");
	  } else {
	    it_optim[nbr_opt] = quadra_breg(beta, xtx, xty, lambda1[m], grd, B);
	  }
	} catch (std::runtime_error& error) {
	  if (verbose > 0) {
	    Rprintf("\nWarning: experiencing numerical instability... ");
	  }
	  if (bullet) {
	    if (verbose > 0) {
	      Rprintf("\nEntering 'bulletproof' mode: switching to proximal algorithm (slower but safer).");
	    } 
	    it_active[m] = 0; // start this lambda all the way back
	    eps = 1e-2 ; // with the proximal settings
	    fun = 1 ;
	    it_optim[nbr_opt] = fista_breg(beta, xtx, xty, grd, lambda1[m], L, pow(eps,2));
	    // reformating the output
	    B    = find(abs(abs(beta) - max(abs(beta))) < ZERO );	    
	  } else {
	    if (verbose > 0) {
	      Rprintf("\nCutting the solution path to this point, as you specified bulletproof=FALSE.");
	    }
	    success_optim = false ;
	  }
	}
      }
      nbr_opt++;
	    
      // _____________________________________________________________
      //
      // (2) OPTIMALITY TESTING
      // _____________________________________________________________

      // dual norm of gradient for unactive variable
      max_grd[m] = as_scalar(sum(abs(grd)) - lambda1[m]) ;
      if (max_grd[m] < 0) {
	max_grd[m] = 0;
      }
      
      // Moving to the next iterate
      it_active[m]++;

      // Cutting the path here if fail to converge or 
      if (!success_optim) {
	break;
      }

      R_CheckUserInterrupt();
    }

    // Record the time ellapsed
    timing[m] = timer.toc() ;

    // Checking convergence status
    if (it_active[m] >= max_iter) {
      converge[m] = 1;
    }
    if (p-B.n_elem > max_feat) {
      converge[m] = 2 ;
    }
    if (!success_optim) {
      converge[m] = 3;
    }

    // Preparing next value of the penalty
    if (converge[m] == 2 || converge[m] == 3) {
      lambda1     =    lambda1.subvec(0,m-1) ;
      converge    =  converge.subvec(0,m)    ;
      max_grd     =   max_grd.subvec(0,m-1)  ;
      it_active   = it_active.subvec(0,m)    ;
      timing      =    timing.subvec(0,m)    ;
      break;
    } else {
      if (naive == 1) {
	coef = join_rows(coef,beta/normx);
      } else {
	coef = join_rows(coef,(1+lambda2)*beta/normx);
      }
      iB = join_cols(iB, m*ones(B.n_elem,1) );
      jB = join_cols(jB, conv_to<mat>::from(B) );
    }
    
  }
  // END OF THE LOOP OVER LAMBDA
  
  return List::create(Named("coefficients") = strans(coef),
		      Named("iB")           = iB          ,
		      Named("jB")           = jB          ,
		      Named("lambda1")      = lambda1     ,
		      Named("it.active")    = it_active   ,
		      Named("it.optim")     = it_optim    ,
		      Named("max.grd")      = max_grd     ,
		      Named("timing")       = timing      ,
		      Named("converge")     = converge    );
}
