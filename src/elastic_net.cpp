/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "elastic_net.h"

using namespace Rcpp;
using namespace arma;

SEXP elastic_net(SEXP X, SEXP XTY, SEXP S, SEXP LAMBDA, SEXP PENSCALE, SEXP XBAR, SEXP YBAR, SEXP NORMX, SEXP NORMY, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP ZERO, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP SPARSE, SEXP MONITOR) {
  
  // Reading input variables
  List   SX       = List       (X)         ; // sparsely encoded SCALED design matrix
  vec    lambda1  = as<vec>    (LAMBDA)    ; // penalty levels
  vec    penscale = as<vec>    (PENSCALE)  ; // penality weights
  vec    xbar     = as<vec>    (XBAR)      ; // mean of the predictors
  double ybar     = as<double> (YBAR)      ; // mean of the response
  vec    normx    = as<vec>    (NORMX)     ; // norm of the predictors
  vec    weights  = as<vec>    (WEIGHTS)   ; // norm of the predictors
  double normy    = as<double> (NORMY)     ; // norm of the predictors
  bool   naive    = as<bool>   (NAIVE)     ; // naive elastic-net or not
  vec    xty      = as<vec>    (XTY)       ; // reponses to predictors vector
  double eps      = as<double> (EPS)       ; // precision required
  double zero     = as<double> (ZERO)      ; // practical zero
  uword  fun      = as<int>    (FUN)       ; // solver (0=quadra, 1=pathwise, 2=fista)
  bool   verbose  = as<bool>   (VERBOSE)   ; // boolean for verbose mode
  bool   sparse   = as<bool>   (SPARSE)    ; // boolean for sparse mode
  int    monitor  = as<int>    (MONITOR)   ; // convergence monitoring (1 == Grandvalet's bound ;-) 2 == Fenchel duality gap)
  uword  max_iter = as<int>    (MAXITER)   ; // max # of iterates of the active set
  uword  max_feat = as<int>    (MAXFEAT)   ; // max # of variables activated
  // Managing the sparse encoding of the structuring matrix
  List SS                                  ; // sparsely encoded structuring matrix
  uvec Si                                  ; // row indices of nonzeros
  uvec Sp                                  ; // col indices of nonzeros
  vec  Sx                                  ; // values of nonzeros
  vec  Sxnzero                             ; // temporary variables used for
  uvec Sjnzero                             ; // updating the Gram matrix
  double lambda2 = 0                       ; // the smooth (ridge) penality
  if (S != R_NilValue) { // Check if their is any structure to read
    SS    = List(S)                        ;
    Si    = as<uvec>(SS[0])                ;
    Sp    = as<uvec>(SS[1])                ;
    Sx    = as<vec> (SS[2])                ;
    lambda2 = Sx[0]                        ;
  }

  // Managing the design matrix is both cases of sparse or dense coding
  mat  x                                   ; // dense encoding of the design matrix
  uvec Xi                                  ; // row indices of nonzeros
  uvec Xp                                  ; // indices of first nonzero of each column
  uvec Xnp                                 ; // # of nonzero in each column
  vec  Xx                                  ; // values of nonzeros  
  uvec j_nz                                ; // consider just the column of X which are non zero
  vec  col_Vx                              ;
  uvec col_Xi                              ;
  vec  col_Xx                              ;
  if (sparse == 1) { // Check how x is encoded for reading
    Xi        = as<uvec>(SX[0]);
    Xp        = as<uvec>(SX[1]); 
    Xnp       = as<uvec>(SX[2]);
    Xx        = as<vec> (SX[3]);
    j_nz = find(Xnp > 0);
  } else { 
    x = as<mat>(SX[0]) - sqrt(weights) * trans(xbar) ;   
  }

  // Initializing "first level" variables (outside of the lambda1 loop)
  uword n        = weights.n_elem        ; // sample size
  uword p        = xty.n_elem            ; // problem size
  uword n_lambda = lambda1.n_elem        ; // # of penalty 
  vec  mu          (n_lambda)            ; // vector of intercept
  uvec A                                 ; // set of currently activated variables
  vec  betaA                             ; // vector of currently activated parameters
  vec  new_col                           ; // column currently added to xtxA
  mat  xtxA                              ; // t(x) * x_A  covariance matrix 
  mat  xAtxA                             ; // t(x_A) * x_A covariance matrix of the activated variable
  vec  xtxw                              ; // t(x_A) * x_A * beta(A)
  mat  R                                 ; // Cholesky decomposition of XAtXA
  vec  grd       = -xty                  ; // smooth part of the gradient
  vec  max_grd   = zeros<vec>(n_lambda)  ; // a vector with the successively reach duality gap
  vec  converge  = zeros<vec>(n_lambda)  ; // a vector indicating if convergence occured (0/1/2)
  uvec it_active = zeros<uvec>(n_lambda) ; // # of loop in the active set for each lambda1
  uvec it_optim                          ; // # of loop in the optimization process for each loop of the active se
  double L0            = 1.0 + lambda2   ; // Lipschitz constant for proximal methods
  vec  timing      (n_lambda)            ; // succesive timing in  
  wall_clock timer                       ; // clock 

  // Initializing "second level" variables (within the active set - for a fixed value of lamdba)
  int         iter     = 0               ; // current iterate
  uword var_in                           ; // currently added variable
  int         nbr_in   = 0               ; // # of currently added variables 
  int         nbr_opt  = 0               ; // # of current calls to the optimization routine
  uvec  are_in   = zeros<uvec>(p)        ; // a vector to check if a variable is already in the active set
  List  out_optim                        ; // the list of output of the optimization function
  uvec  null                             ; // stores the variables which go to zero during optimization
  uword swap                             ; // stores the variable whose sign has swaped during optimization
  vec   grd_norm (p)                     ; // current value of the grd_norm for each variable 
  vec   pen      (p)                     ; // current vector of penalties 
  mat   nonzeros                         ; // contains non-zero value of beta
  mat   iA                               ; // contains row indices of the non-zero values
  mat   jA                               ; // contains column indices of the non-zero values

  // Additional variable for convergence monitoring
  vec delta_hat    ;
  vec delta_star   ;
  vec gamma        ;
  uvec Ac          ;
  double nu        ; 
  double quad_loss ;
  vec J_hat        ;
  mat J_star       ;

  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (int m=0; m<n_lambda; m++) {
    if (verbose == 1) {Rprintf("\n lambda1 = %f",lambda1(m)) ;}
    // _____________________________________________________________
    //
    // START THE ACTIVE SET ALGORITHM
    // _____________________________________________________________
    //
    
    // current vector of weighted penalities
    pen = lambda1[m]*penscale ;
    // dual norm of gradient for unactive variable
    grd_norm = abs(grd) - pen ;
    // dual norm of gradient for active variables
    grd_norm.elem(A) = abs(grd.elem(A) + pen.elem(A) % signs(betaA, zero)) ;
    // variable associated with the highest optimality violation
    max_grd[m] = grd_norm.max(var_in) ;
    if (max_grd[m] < 0) {
      max_grd[m] = 0;
    }
     
    while (max_grd[m] > eps & it_active[m] < max_iter) {
      // _____________________________________________________________
      //
      // (1) VARIABLE ACTIVATION IF APPLICABLE
      // _____________________________________________________________

      // Check if the variable is already in the active set
      if (are_in[var_in] == 0) {
	// If this is a newly added variable, then
	A.resize(nbr_in+1)     ; // update the active set
	A[nbr_in] = var_in     ;
	betaA.resize(nbr_in+1) ; // update the vector of active parameters
	betaA[nbr_in]  = 0.0   ;
	
	// Update the xtxA and xAtxA matrices
	if (nbr_in > 0) {
	  xAtxA = join_cols(xAtxA, xtxA.row(var_in)) ;
	}

	if (sparse == 1) {
	  // if any nonzero in the X[, var] column, do the sparse product
	  new_col = zeros<vec>(p) ;
	  if (Xnp[var_in] > 0) {
	    col_Vx = zeros<vec>(n) ;
	    col_Vx.elem(Xi.subvec(Xp[var_in],Xp[var_in+1]-1)) = Xx.subvec(Xp[var_in],Xp[var_in+1]-1);
    	    // loop along each column of X
	    for (int j=0; j<j_nz.n_elem; j++) {
	      col_Xx = Xx.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
	      col_Xi = Xi.subvec(Xp[j_nz[j]],Xp[j_nz[j]+1]-1) ;
	      new_col[j_nz[j]] = dot(col_Xx, col_Vx.elem(col_Xi));
	    }
	  }
	  new_col = new_col - n * xbar * as_scalar(xbar[var_in]);
	} else {
	  new_col = trans(x) * x.col(var_in);
	}
	
	if (lambda2 > 0) { 
	  // Sparse product with the structurating matrix
	  Sxnzero = Sx.subvec(Sp[var_in],Sp[var_in+1]-1) ;
	  Sjnzero = Si.subvec(Sp[var_in],Sp[var_in+1]-1) ;
	  new_col.elem(Sjnzero) +=  new_col.elem(Sjnzero) % Sxnzero ;
	}
	
	xtxA  = join_rows(xtxA, new_col) ;
	xAtxA = join_rows(xAtxA, trans(xtxA.row(var_in))) ;
	if (fun == 0) {
	  R = cholupdate(R, xAtxA) ;
	}
	if (fun == 1) {
	  xtxw.resize(nbr_in+1) ;
	  xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
	}
	if (verbose == 1) {Rprintf("newly added variable %i\n",var_in);}
	are_in[var_in] = 1;
	nbr_in++;
      } else {
	if (verbose == 1) {Rprintf("already in %i\n",var_in);}
      }
      
      // _____________________________________________________________
      //
      // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      // _____________________________________________________________
      //      
      switch (fun) {
      case 1 :
	out_optim = pathwise_enet(betaA, xAtxA, xty.elem(A), xtxw, pen.elem(A), lambda2, eps);
	xtxw = as<vec>(out_optim[2]);
	break;
      case 2 :
	out_optim = fista_lasso(betaA, xAtxA, xty.elem(A), pen.elem(A), L0, eps);
	L0 = as<double>(out_optim[2]);
	break;
      default:
	out_optim = quadra_enet(betaA, R, xty.elem(A), signs(grd.elem(A), zero), pen.elem(A), zero);
      }
      // save the number of iterates performed along the optimization process
      it_optim.reshape(nbr_opt + 1,1) ;
      it_optim[nbr_opt] = as<int>(out_optim[0]) + 1 ;
      nbr_opt++;

      // _____________________________________________________________
      //
      // (3) VARIABLE DELETION IF APPLICABLE
      // _____________________________________________________________
      //
      if (fun == 0) {
	if (out_optim.size() > 2) {
	  betaA = as<vec>(out_optim[2]);
	  swap  = as<uword>(out_optim[3]);
	  grd = -xty + xtxA * betaA;
	  if (fabs(grd(A(swap)) + pen(A(swap)) * sign(betaA(swap), zero)) > eps) {
	    betaA  = as<vec>(out_optim[1]);
	    grd = -xty + xtxA * betaA;
	    if (verbose == 1) {Rprintf("removing variable %i\n",swap);}
	    are_in[A(swap)]  = 0 ;
	    A.shed_row(swap)     ;
	    betaA.shed_row(swap) ;
	    xtxA.shed_col(swap)  ;
	    xAtxA.shed_col(swap) ;
	    xAtxA.shed_row(swap) ;
	    R = choldowndate(R, swap) ;
	    nbr_in--;
	  }
	} else {
	  betaA  = as<vec>(out_optim[1]) ;
	  grd = -xty + xtxA * betaA;
	  // check sign cohenrency between gradient and parameter lastly added
	  if (sign(grd[var_in], zero) == sign(betaA[nbr_in-1], zero)) {
	    out_optim = quadra_enet(betaA, R, xty.elem(A), signs(grd.elem(A), zero), pen.elem(A), zero);
	    betaA  = as<vec>(out_optim[1]) ;
	    grd = -xty + xtxA * betaA;
	    // add the corresponding iteration
	    it_optim.reshape(nbr_opt + 1,1) ;
	    it_optim[nbr_opt] = it_optim[nbr_opt] + out_optim[0] + 1 ;
	  }
	}
      } else {
	betaA  = as<vec>(out_optim[1]) ;
	grd = -xty + xtxA * betaA;
	null = sort(find(abs(betaA) + (abs(grd.elem(A)) - pen.elem(A)) < zero),1) ;
	if (!null.is_empty()) {
	  for (int j=0; j<null.n_elem; j++) {
	    if (verbose == 1) {Rprintf("removing variable %i\n",null[j]);}
	    are_in[null[j]] = 0     ;
	    A.shed_row(null[j])     ;
	    betaA.shed_row(null[j]) ;
	    xtxA.shed_col(null[j])  ;
	    xAtxA.shed_col(null[j]) ;
	    xAtxA.shed_row(null[j]) ;
	    if (fun == 0) {
	      R = choldowndate(R, null[j]);
	    }
	    if (fun == 1) {
	      xtxw.shed_row(null[j]);
	    }
	    nbr_in--;
	  }
	}
      }
      
      // _____________________________________________________________
      //
      // (4) OPTIMALITY TESTING
      // _____________________________________________________________
      
      // dual norm of gradient for unactive variable
      grd_norm = abs(grd) - pen ;
      // dual norm of gradient for active variables
      grd_norm.elem(A) = abs(grd.elem(A) + pen.elem(A) % signs(betaA, zero)) ;
      // variable associated with the highest optimality violation
      max_grd[m]  = grd_norm.max(var_in) ;
      if (max_grd[m] < 0) {
	max_grd[m] = 0;
      }
      
      if (monitor > 0) {
	// _____________________________________________________________
	//
	// FOLLOWING CONVERGENCE BY COMPLETE MONITORING
	// _____________________________________________________________

	// gamma equals the max |gradient|
	gamma = pen % grd/lambda1[m] ;
	nu = norm(gamma, "inf");
	
 	J_hat.reshape(nbr_opt,1) ;
	quad_loss =  pow(normy,2) + dot(betaA,xAtxA * betaA) - 2*dot(betaA, xty.elem(A)) ;
	J_hat[nbr_opt-1] =  0.5*quad_loss - dot(betaA, gamma.elem(A));
	
	delta_hat.reshape(nbr_opt,1) ;
	if (monitor == 1) {
	  Ac = find(gamma > nu); // set of adversarial variables outside the boundary
	  // Grandvalet's bound
	  delta_hat[nbr_opt-1] = J_hat[nbr_opt-1] - (lambda1[m]/nu) * J_hat[nbr_opt-1]  - (pow(lambda1[m],2)/(2*lambda2))*((lambda1[m]*(p-Ac.n_elem))/nu + pow(norm(gamma.elem(Ac),2)/nu,2) - p);
	} else {
	  // Fenchel's bound
	  if (nu < lambda1[m]) {
	    nu <- lambda1[m];
	  }
	  delta_hat[nbr_opt-1] = 0.5 * quad_loss * (1+pow(lambda1[m]/nu,2)) + sum(abs(pen.elem(A) % betaA)) + (lambda1[m]/nu)*(dot(betaA,xty.elem(A))-pow(normy,2));
	}
	// keep the smallest bound reached so far for a given lambda value
 	if (nbr_opt>1) {
	  if (J_hat[nbr_opt-2] < J_hat[nbr_opt-1]) {
	    if (delta_hat[nbr_opt-1] > delta_hat[nbr_opt-2] - (J_hat[nbr_opt-2] - J_hat[nbr_opt-1])) {
	      delta_hat[nbr_opt-1] = delta_hat[nbr_opt-2];
	    }
	  }
   	}
      }
      
      // Moving to the next iterate
      it_active[m]++;
      
      R_CheckUserInterrupt();    
    }
    
    // the reference parameter (obtained once optimum is met)
    if (monitor > 0) {
      if (it_active[m] > 0) {
	J_star = join_cols(J_star, ones(it_active[m],1) * J_hat[nbr_opt-1]) ;
      }
    }
    
    // Preparing next value of the penalty
    if (naive == 1) {
      nonzeros = join_cols(nonzeros, betaA/normx.elem(A));
    } else {
      nonzeros = join_cols(nonzeros, (1+lambda2)*betaA/normx.elem(A));
    }
    iA = join_cols(iA, m*ones(betaA.n_elem,1) );
    jA = join_cols(jA, conv_to<mat>::from(A) ) ;
    
    mu[m] = ybar - dot(xbar.elem(A), (1+lambda2)*betaA);
    
    if (it_active[m] >= max_iter) {
      converge[m] = 1;
    }

    timing[m] = timer.toc() ;
    if (nbr_in > max_feat) {
      converge[m] = 2 ;
      lambda1      =    lambda1.subvec(0,m) ;
      converge    =  converge.subvec(0,m) ;
      mu          =        mu.subvec(0,m) ;
      max_grd     =   max_grd.subvec(0,m) ;
      it_active   = it_active.subvec(0,m) ;
      timing      =    timing.subvec(0,m) ;
      break;
    }
    
  }
  // END OF THE LOOP OVER LAMBDA
  
  // Updating monitored quantities
  if (monitor > 0) {
    delta_star = J_hat - J_star;
  }
  
  return List::create(Named("mu")         = mu        ,
			    Named("nzeros")     = nonzeros  ,
			    Named("iA")         = iA        ,
			    Named("jA")         = jA        ,
 			    Named("lambda1")    = lambda1   ,
 			    Named("nbr.in")     = nbr_in    ,
 			    Named("it.active")  = it_active ,
 			    Named("it.optim")   = it_optim  ,
 			    Named("max.grd")    = max_grd   ,
 			    Named("timing")     = timing    ,
 			    Named("delta.hat")  = delta_hat ,
 			    Named("delta.star") = delta_star,
 			    Named("converge")   = converge  );
    
}
