/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */

#include "elastic_net.h"

using namespace Rcpp;
using namespace arma;

SEXP elastic_net(SEXP X, SEXP XTY, SEXP S, SEXP LAMBDA1, SEXP LAMBDA2, SEXP XBAR, SEXP NORMX, SEXP NORMY, SEXP WEIGHTS, SEXP NAIVE, SEXP EPS, SEXP MAXITER, SEXP MAXFEAT, SEXP FUN, SEXP VERBOSE, SEXP SPARSE, SEXP MONITOR) {

  // disable messages being printed to the err2 stream
  std::ostream nullstream(0);
  set_stream_err2(nullstream);

  // Reading input variables
  List   SX       = List       (X)         ; // sparsely encoded SCALED design matrix
  vec    lambda1  = as<vec>    (LAMBDA1)   ; // penalty levels
  double lambda2  = as<double> (LAMBDA2)   ; // the smooth (ridge) penality
  vec    xbar     = as<vec>    (XBAR)      ; // mean of the predictors
  vec    normx    = as<vec>    (NORMX)     ; // norm of the predictors
  vec    weights  = as<vec>    (WEIGHTS)   ; // norm of the predictors
  double normy    = as<double> (NORMY)     ; // norm of the predictors
  bool   naive    = as<bool>   (NAIVE)     ; // naive elastic-net or not
  vec    xty      = as<vec>    (XTY)       ; // reponses to predictors vector
  double eps      = as<double> (EPS)       ; // precision required
  uword  fun      = as<int>    (FUN)       ; // solver (0=quadra, 1=pathwise, 2=fista)
  int    verbose  = as<int>    (VERBOSE)   ; // int for verbose mode (0/1/2)
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
  if (S != R_NilValue) { // Check if their is any structure to read
    SS    = List(S)                        ;
    Si    = as<uvec>(SS[0])                ;
    Sp    = as<uvec>(SS[1])                ;
    Sx    = as<vec> (SS[2])                ;
  }

  // Managing the design matrix is both cases of sparse or dense coding
  mat  x                                   ; // dense encoding of the design matrix
  uvec Xi                                  ; // row indices of nonzeros
  uvec Xp                                  ; // indices of first nonzero of each column
  uvec Xnp                                 ; // # of nonzero in each column
  vec  Xx                                  ; // values of nonzeros
  uvec j_nz                                ; // consider just the column of X which are non zero
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
  uvec A                                 ; // set of currently activated variables
  vec  betaA                             ; // vector of currently activated parameters
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
  uword var_in                           ; // currently added variable
  int   nbr_in   = 0                     ; // # of currently added variables
  int   nbr_opt  = 0                     ; // # of current calls to the optimization routine
  uvec  are_in   = zeros<uvec>(p)        ; // a vector to check if a variable is already in the active set
  List  out_optim                        ; // the list of output of the optimization function
  bool  success_optim = true             ; // was the internal system resolution successful?
  uvec  null                             ; // stores the variables which go to zero during optimization
  vec   grd_norm (p)                     ; // current value of the grd_norm for each variable
  mat   nonzeros                         ; // contains non-zero value of beta
  mat   iA                               ; // contains row indices of the non-zero values
  mat   jA                               ; // contains column indices of the non-zero values

  // Additional variable for convergence monitoring
  vec D_hat    ;
  vec D_star   ;
  vec J_hat    ;
  mat J_star   ;

  // _____________________________________________________________
  //
  // START THE LOOP OVER LAMBDA
  timer.tic();
  for (int m=0; m<n_lambda; m++) {
    if (verbose == 2) {Rprintf("\n lambda1 = %f",lambda1(m)) ;}
    // _____________________________________________________________
    //
    // START THE ACTIVE SET ALGORITHM
    // _____________________________________________________________
    //

    // dual norm of gradient for unactive variable
    grd_norm = abs(grd) - lambda1[m] ;
    // dual norm of gradient for active variables
    grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * signs(betaA)) ;
    // variable associated with the highest optimality violation
    max_grd[m] = grd_norm.max(var_in) ;
    if (max_grd[m] < 0) {max_grd[m] = 0;}

    while ((max_grd[m] > eps) && (it_active[m] < max_iter)) {
      // _____________________________________________________________
      //
      // (1) VARIABLE ACTIVATION IF APPLICABLE
      // _____________________________________________________________

      // Check if the variable is already in the active set
      if (are_in[var_in] == 0) {
	add_var_enet(n, nbr_in, var_in, betaA, A, x, xtxA, xAtxA, xtxw, R, lambda2, xbar, Xx, Xi, Xp, Xnp, j_nz, Sx, Si, Sp, sparse, fun) ;
	if (verbose == 2) {Rprintf("newly added variable %i\n",var_in);}
	are_in[var_in] = 1;
	nbr_in++;
      } else {
	if (verbose == 2) {Rprintf("already in %i\n",var_in);}
      }

      // _____________________________________________________________
      //
      // (2) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
      // _____________________________________________________________
      //

      it_optim.reshape(nbr_opt + 1,1) ;
      switch (fun) {
      case 1 :
	it_optim[nbr_opt] = pathwise_enet(betaA, xAtxA, xty.elem(A), xtxw, lambda1[m], null, lambda2, pow(eps,2));
	break;
      case 2 :
	it_optim[nbr_opt] = fista_lasso(betaA, xAtxA, xty.elem(A), lambda1[m], null, L0, pow(eps,2));
	break;
      default:
	try {
	  it_optim[nbr_opt] = quadra_enet(betaA, R, xty.elem(A), signs(grd.elem(A)), lambda1[m], null);
	} catch (std::runtime_error &error) {
	  if (verbose > 0) {
	    Rprintf("\nWarning: singular system at this stage of the solution path, cutting here.\n");
	  }
	  success_optim = false ;
	}
      }
      // update the smooth part of the gradient
      grd = -xty + xtxA * betaA;
      nbr_opt++;

      // _____________________________________________________________
      //
      // (3) VARIABLE DELETION IF APPLICABLE
      // _____________________________________________________________
      //
      // removing variables zeroed during optimization
      if (!null.is_empty()) {
	if (verbose == 2) {
	  for (int j=0; j<null.n_elem; j++) {Rprintf("removing variable %i\n",null[j]);}
	}
	remove_var_enet(nbr_in,are_in,betaA,A,xtxA,xAtxA,xtxw,R,null,fun) ;
      }

      // Seems useless (handled by the while loop over the optimality conditions)
      // For Quadratic solver only: check sign coherency between the updated gradient and variable lastly optimized
      // if (fun == 0 & it_optim[nbr_opt-1] <= 1) {
      // 	// if not, reoptimize...
      // 	if (nbr_in > 0 & sign(grd[var_in], ZERO) == sign(betaA[nbr_in-1], ZERO)) {
      // 	  it_optim.reshape(nbr_opt + 1,1) ;
      // 	  it_optim[nbr_opt] += quadra_enet(betaA, R, xty.elem(A), signs(grd.elem(A), ZERO), lambda1[m], null, ZERO);
      // 	  grd = -xty + xtxA * betaA;
      // 	}
      // }

      // _____________________________________________________________
      //
      // (4) OPTIMALITY TESTING
      // _____________________________________________________________

      // dual norm of gradient for unactive variable
      grd_norm = abs(grd) - lambda1[m] ;
      // dual norm of gradient for active variables
      grd_norm.elem(A) = abs(grd.elem(A) + lambda1[m] * signs(betaA)) ;
      // variable associated with the highest optimality violation
      max_grd[m]  = grd_norm.max(var_in) ;
      if (max_grd[m] < 0) {max_grd[m] = 0;}

      if (monitor > 0) {
	// _____________________________________________________________
	//
	// (OPTIONAL) FOLLOWING CONVERGENCE BY COMPLETE MONITORING
	// _____________________________________________________________
	bound_to_optimal(betaA, xAtxA, xty, grd, lambda1[m], lambda2, normy, A, monitor, J_hat, D_hat) ;
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

    // Record the time ellapsed
    timing[m] = timer.toc() ;

    // Checking convergence status
    if (it_active[m] >= max_iter) {
      converge[m] = 1;
    }
    if (nbr_in > max_feat) {
      converge[m] = 2 ;
    }
    if (!success_optim) {
      converge[m] = 3;
    }

    // Stop now if relevant
    if (converge[m] == 2 || converge[m] == 3) {
      lambda1     =    lambda1.subvec(0,m-1) ;
      converge    =  converge.subvec(0,m)    ;
      max_grd     =   max_grd.subvec(0,m-1)  ;
      it_active   = it_active.subvec(0,m)    ;
      timing      =    timing.subvec(0,m)    ;
      break;
    } else {
      // Preparing next value of the penalty
      if (naive == 1) {
	nonzeros = join_cols(nonzeros, betaA/normx.elem(A));
      } else {
	nonzeros = join_cols(nonzeros, (1+lambda2)*betaA/normx.elem(A));
      }
      iA = join_cols(iA, m*ones(betaA.n_elem,1) );
      jA = join_cols(jA, conv_to<mat>::from(A) ) ;
    }
  }
  // END OF THE LOOP OVER LAMBDA

  // Updating monitored quantities
  if (monitor > 0) {
    D_star = J_hat - J_star;
  }

  return List::create(Named("nzeros")     = nonzeros ,
		      Named("iA")         = iA       ,
		      Named("jA")         = jA       ,
		      Named("lambda1")    = lambda1  ,
		      Named("nbr.in")     = nbr_in   ,
		      Named("it.active")  = it_active,
		      Named("it.optim")   = it_optim ,
		      Named("max.grd")    = max_grd  ,
		      Named("timing")     = timing   ,
		      Named("delta.hat")  = D_hat    ,
		      Named("delta.star") = D_star   ,
		      Named("converge")   = converge );

}

void cholupdate(mat &R , mat &XtX) {
  int p = XtX.n_cols;

  if (p == 1) {
    R = sqrt(XtX);
  } else {
    colvec rp  = zeros<colvec>(p,1);
    rp.subvec(0,p-2) = solve (trimatl(strans(R)), XtX.submat(0,p-1,p-2,p-1));
    rp(p-1) = sqrt(XtX(p-1,p-1) - dot(rp,rp));
    R = join_rows( join_cols(R, zeros<mat>(1,p-1)) , rp);
  }
}

void choldowndate(mat &R, int j) {

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
}

void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, mat &x, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, double &lambda2,
		  vec &xbar, vec &Xx, uvec &Xi, uvec &Xp, uvec &Xnp, uvec &j_nz, vec &Sx, uvec &Si, uvec& Sp, bool &sparse, uword &fun) {

  vec  Sxnzero   ; // temporary variables used for
  uvec Sjnzero   ; // updating the Gram matrix
  vec  new_col   ; // column currently added to xtxA
  vec  col_Vx    ;
  uvec col_Xi    ;
  vec  col_Xx    ;

  // If this is a newly added variable, then
  A.resize(nbr_in+1)     ; // update the active set
  A[nbr_in] = var_in     ;
  betaA.resize(nbr_in+1) ; // update the vector of active parameters
  betaA[nbr_in]  = 0.0   ;

  // Update the xtxA and xAtxA matrices
  if (nbr_in > 0) {
    xAtxA = join_cols(xAtxA, xtxA.row(var_in)) ;
  }

  // This would greatly simply if Iused the armadillo sparse features
  // Yet, it requires importation of sparseMatrix fro R vuia Rcpp
  if (sparse == 1) {
    // if any nonzero in the X[, var] column, do the sparse product
    new_col = zeros<vec>(xbar.n_elem) ;
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
    new_col.elem(Sjnzero) += new_col.elem(Sjnzero) % Sxnzero ;
  }

  xtxA  = join_rows(xtxA, new_col) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.row(var_in))) ;
  if (fun == 0) {
    cholupdate(R, xAtxA) ;
  }
  if (fun == 1) {
    xtxw.resize(nbr_in+1) ;
    xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
  }

}

void remove_var_enet(int &nbr_in, uvec &are_in, vec &betaA, uvec &A, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, uvec &null, uword &fun) {

  for (int j=0; j<null.n_elem; j++) {
    are_in[A(null[j])]  = 0 ;
    A.shed_row(null[j])     ;
    betaA.shed_row(null[j]) ;
    xtxA.shed_col(null[j])  ;
    xAtxA.shed_col(null[j]) ;
    xAtxA.shed_row(null[j]) ;
    if (fun == 1) {
      xtxw.shed_row(null[j]);
    }
    if (fun == 0) {
      choldowndate(R, null[j]) ;
    }
    nbr_in--;
  }

}

void bound_to_optimal(vec &betaA,
		      mat &xAtxA,
		      vec &xty,
		      vec &grd,
		      double &lambda1,
		      double &lambda2,
		      double &normy,
		      uvec &A,
		      int &monitor,
		      vec &J_hat,
		      vec &D_hat) {

  // to store the results
  int dim = J_hat.n_elem ;
  J_hat.resize(dim+1);
  D_hat.resize(dim+1);

  // gamma equals the max |gradient|
  vec gamma = grd ;
  double nu = norm(gamma, "inf");
  int p = xty.n_elem ;

  double quad_loss =  pow(normy,2) + dot(betaA,xAtxA * betaA) - 2*dot(betaA, xty.elem(A)) ;
  J_hat(dim)  =  0.5*quad_loss - dot(betaA, gamma.elem(A));

  if (monitor == 1) {
    uvec Ac = find(gamma > nu); // set of adversarial variables outside the boundary
    // Grandvalet's bound
    D_hat(dim) = J_hat(dim) - (lambda1/nu) * J_hat(dim) - (pow(lambda1,2)/(2*lambda2))*((lambda1*(p-Ac.n_elem))/nu + pow(norm(gamma.elem(Ac),2)/nu,2)-p);
  } else {
    // Fenchel's bound
    if (nu < lambda1) {
      nu = lambda1;
    }
    D_hat(dim) = 0.5 * quad_loss * (1+pow(lambda1/nu,2)) + sum(abs(lambda1*betaA)) + (lambda1/nu)*(dot(betaA,xty.elem(A))-pow(normy,2));
  }

  // keep the smallest bound reached so far for a given lambda value
  if (dim>0) {
    if (J_hat[dim-1] < J_hat[dim]) {
      if (D_hat[dim] > D_hat[dim-1] - (J_hat[dim-1] - J_hat[dim])) {
	D_hat[dim] = D_hat[dim-1];
      }
    }
  }

}
