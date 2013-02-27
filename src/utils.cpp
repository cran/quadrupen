#include "utils.h"

using namespace Rcpp;
using namespace arma;

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


void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, mat &x, mat &xt, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, double &lambda2, vec &xbar, sp_mat &spS, bool &usechol, uword &fun) {

  vec  new_col   ; // column currently added to xtxA

  A.resize(nbr_in+1)     ; // update the active set
  A[nbr_in] = var_in     ;
  betaA.resize(nbr_in+1) ; // update the vector of active parameters
  betaA[nbr_in]  = 0.0   ;

  new_col = xt * x.col(var_in);
  if (lambda2 > 0) {
    // Adding the column corresponding to the structurating matrix
    new_col += spS.col(var_in);
  }

  // UPDATE THE xtxA AND xAtxA MATRICES
  if (nbr_in > 0) {
    xAtxA = join_cols(xAtxA, xtxA.row(var_in)) ;
  }
  xtxA  = join_rows(xtxA, new_col) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.row(var_in))) ;

  if (fun == 0 & usechol == 1) {
    cholupdate(R, xAtxA) ;
  }

  if (fun == 1) {
    xtxw.resize(nbr_in+1) ;
    xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
  }
}

void add_var_enet(uword &n, int &nbr_in, uword &var_in, vec &betaA, uvec &A, sp_mat &x, sp_mat &xt, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, double &lambda2, vec &xbar, sp_mat &spS, bool &usechol, uword &fun) {

  vec  new_col   ; // column currently added to xtxA

  A.resize(nbr_in+1)     ; // update the active set
  A[nbr_in] = var_in     ;
  betaA.resize(nbr_in+1) ; // update the vector of active parameters
  betaA[nbr_in]  = 0.0   ;

  new_col = xt * x.col(var_in) - n * xbar * as_scalar(xbar[var_in]);
  if (lambda2 > 0) {
    // Adding the column corresponding to the structurating matrix
    new_col += spS.col(var_in);
  }

  // UPDATE THE xtxA AND xAtxA MATRICES
  if (nbr_in > 0) {
    xAtxA = join_cols(xAtxA, xtxA.row(var_in)) ;
  }
  xtxA  = join_rows(xtxA, new_col) ;
  xAtxA = join_rows(xAtxA, trans(xtxA.row(var_in))) ;

  if (fun == 0 & usechol == 1) {
    cholupdate(R, xAtxA) ;
  }

  if (fun == 1) {
    xtxw.resize(nbr_in+1) ;
    xtxw(nbr_in) = dot(xAtxA.col(nbr_in),betaA);
  }
}

void remove_var_enet(int &nbr_in, uvec &are_in, vec &betaA, uvec &A, mat &xtxA, mat &xAtxA, mat &xtxw, mat &R, uvec &null, bool &usechol, uword &fun) {

  for (int j=0; j<null.n_elem; j++) {
    are_in[A(null[j])]  = 0 ;
    A.shed_row(null[j])     ;
    betaA.shed_row(null[j]) ;
    if (fun == 1) {
      xtxw.shed_row(null[j]);
    }
    xtxA.shed_col(null[j])  ;
    xAtxA.shed_col(null[j]) ;
    xAtxA.shed_row(null[j]) ;
    if (fun == 0 & usechol == 1) {
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

vec cg(mat A, vec b, vec x, double tol) {

  vec r = b - A * x;
  vec p = r ;
  double rs_old = sum(square(r)) ;

  double rs_new = rs_old ;
  int i = 0;
  double alpha ;
  mat Ap ;

  while (sqrt(rs_new) > tol & i < 1e3) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    x += alpha * p ;
    r -= alpha * Ap ;
    // Polak-Ribière update
    rs_new = dot(r,-alpha * Ap);
    p = r + rs_new/rs_old*p;
    rs_old = rs_new;
    i++;
  }

  // Rprintf("\n nb of iterate %d",i);
  return(x);
}

// Can't find a reasonable Preconditioner that does not
// require a computational burden equivalent to a Cholesky decomposition
vec pcg(mat A, mat P, vec b, vec x, double tol) {

  vec r = b - A * x;
  vec z = P * r;
  vec p = z ;
  //double rs_old = sum(square(r)) ;
  double rs_old = dot(r,z) ;

  double rs_new = rs_old ;
  int i = 0;
  double alpha ;
  mat Ap ;

  while (sqrt(rs_new) > tol & i < 1e3) {
    Ap = A * p;
    alpha = rs_old/dot(p,Ap) ;
    x += alpha * p ;
    r -= alpha * Ap ;
    // Polak-Ribière update
    z = P * r ;
    rs_new = dot(z,-alpha * Ap);
    p = z + rs_new/rs_old*p;
    rs_old = rs_new;
    i++;
  }

  // Rprintf("\n nb of iterate %d",i);
  return(x);
}
