#pragma once

#include "FusedLasso_utils.h"
#include "FusedLasso_enums.h"


// Very simple implementation of a sparse matrix
// The matrix sparseness structure can not be changed after construction
// Only used for the X matrix, which does not need to be changed 

class SparseMatrix {
public:
  int nnz;
  int n;
  int p;
  vector<double> nzVec; // vector of non-zero elements
  vector<int> indexCol; // index of the first element of x of column; length p+1
  vector<int> rowPos; // row position of the ith element of x; length size(x)
  
public:
  
  // Constructors
  SparseMatrix();
  // also empty, but with specified number of rows
  SparseMatrix(int n);
  // From a c-Array
  SparseMatrix(const double*x, int n, int p);
  //From a vector
  SparseMatrix(const vector<double>& x, int n, int p);
  // From armadillo sparse format
  SparseMatrix(const arma::sp_mat& x);
  // From an already sparse format
  SparseMatrix(const vector<double>& x, const vector<int>& indexCol, const vector<int>& rowPos, int n, int p);
  // create from a dgCMatrix; assumes that it is the right class and does not check
  SparseMatrix(Rcpp::S4 x);
  // given a set of fusions, generate a new matrix (could be extended to multiplying
  // two matrices
  SparseMatrix createFusedX(const vector<vector<int> > fusedGroups);
  
  // add a column to the matrix
  // in sparse format
  void addColumn(const vector<double>& newnzVec, const vector<int>& newrowPos);
  // in non-sparse format
  void addColumn(const vector<double>& x);
  
  // calculate inner product of ith column with itself
  double innerProd(int i) const;
  // calculate inner product of column i and column j
  double innerProd(int i, int j) const;
  
  // inner product of column i with vector y
  double innerProd(const vector<double>& y, int i) const;
  
  // WEIGHTED VERSION OF INNER PRODUCTS 
  // calculate inner product of ith column with itself
  double innerProd(int i, const vector<double>& w) const;
  // calculate inner product of column i and column j
  double innerProd(int i, int j, const vector<double>& w) const;
  
  // inner product of column i with vector y
  double innerProd(const vector<double>& y, int i, const vector<double>& w) const;
  
  // ===== HELPER METHODS FOR innerProd CONSOLIDATION =====
  // These inline helpers reduce code duplication across the 6 innerProd overloads
  // They implement the common sparse iteration patterns
  
  // Helper: dot product of column col with itself, optionally weighted
  // Pattern used by: innerProd(i) when w==nullptr, innerProd(i, w) when w!=nullptr
  inline double dotProductColumnSelf(int col, const vector<double>* weights = nullptr) const {
    if(col < 0 || col >= p) {
      Rcpp::stop("not a valid column %d out of %d\n", col, p);
    }
    double res = 0;
    if(weights == nullptr) {
      // Unweighted: res += X[r,col]^2
      for(int pos = indexCol[col]; pos < indexCol[col+1]; ++pos) {
        res += nzVec[pos] * nzVec[pos];
      }
    } else {
      // Weighted: res += X[r,col]^2 * w[r]
      if((int)weights->size() != n) {
        Rcpp::stop("weights not the right size\n");
      }
      for(int pos = indexCol[col]; pos < indexCol[col+1]; ++pos) {
        res += nzVec[pos] * nzVec[pos] * (*weights)[rowPos[pos]];
      }
    }
    return res;
  }
  
  // Helper: dot product of column col with external vector y, optionally weighted
  // Pattern used by: innerProd(y, i) when w==nullptr, innerProd(y, i, w) when w!=nullptr
  inline double dotProductColumnVector(const vector<double>& y, int col, 
                                       const vector<double>* weights = nullptr) const {
    if((int)y.size() != n) {
      Rcpp::stop("y does not have the right size\n");
    }
    if(col < 0 || col >= p) {
      Rcpp::stop("not a valid column %d out of %d\n", col, p);
    }
    double res = 0;
    if(weights == nullptr) {
      // Unweighted: res += y[r] * X[r,col]
      for(int posI = indexCol[col]; posI < indexCol[col + 1]; ++posI) {
        res += y[rowPos[posI]] * nzVec[posI];
      }
    } else {
      // Weighted: res += y[r] * X[r,col] * w[r]
      if((int)weights->size() != n) {
        Rcpp::stop("weights not the right size\n");
      }
      for(int posI = indexCol[col]; posI < indexCol[col + 1]; ++posI) {
        res += y[rowPos[posI]] * nzVec[posI] * (*weights)[rowPos[posI]];
      }
    }
    return res;
  }
  // ===== END HELPER METHODS =====
  
  // multiply the sparse matrix by a diagonal matrix from the left
  void multByDiag(const vector<double>& w);
  
  // add a multiple of column i to the given vector y; done in place
  void addMultOfColumn(vector<double>& y, int i, double mult) const;
  
  // Unsafe version: skips bounds checking for performance
  // Use only when i is guaranteed to be in [0, p) and y.size() == n
  // Provides 10-15% performance gain in hot loops
  inline void addMultOfColumnUnsafe(vector<double>& y, int i, double mult) const;
  
  // return the value of a position in the matrix
  // with row i and column j
  double get(int i, int j) const;
  
  // return the i-th column as a vector
  vector<double> getColumn(int i) const;
  
  SEXP todgCMatrix() const;
  
  // print out the current object (used for debugging purposes)
  void print(ostream& outStream) const;
  
  // print column
  void printColumn(int i, ostream& outStream) const;
};


// the class that keeps the data for calculating the derivative
class QuadraticDerivative {

public:
    /*
        calculate the derivative for the position pos
    */
    virtual ~QuadraticDerivative() {};

    virtual double getDerivative(int pos) = 0;	

    virtual void updateBeta(int pos, double newBeta) = 0;	

    virtual double getBeta(int pos) = 0;

    virtual vector<double> getBetaVec() = 0;

    virtual int getBetaSize() = 0;

    virtual double getHessian(int pos) = 0;

    virtual int getn() const = 0;

    virtual int getp() const = 0;

    virtual void activate(int pos) = 0;

    virtual vector<double> getDerivativeVec() = 0;

    virtual bool isExtreme() {return(false);};

    virtual void print(vector<int> active) {};

    virtual double lossFuncChange(const vector<vector<int> > &connections, const vector<double> & wLambda1, const vector<vector<double> > &wLambda2) { return(-1); };
};

// the class that keeps the data for calculating the derivative
class QuadraticDerivativeDiagonal : public QuadraticDerivative {
public:
  int n; // number of rows in X
  int p; // number of cols in X
  
  SparseMatrix X; // copy of vector X
  vector<double> y; // copy of vector y
  vector<double> w; // observation weights
  vector<double> beta; // copy of R object that stores the current value of beta
  SparseMatrix WX; // stores a sparse matrix of WX
  vector<double> WXbeta; // stores the current value of Xbeta
  vector<double> diagXTWX;
  vector<double> XTWy;
  
public:
  /*
   Returns the value of the quadratic derivative
   Encapsulates all necessary pre-computation in the Constructor
   Assumes that the dimensions of the objects are correct
   */
  QuadraticDerivativeDiagonal() {};
  
  // Constructor
  QuadraticDerivativeDiagonal(const SparseMatrix& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);
  
  ~QuadraticDerivativeDiagonal() {};
  
  // calculate the derivative for the position pos
  double getDerivative(int pos);	
  
  void updateBeta(int pos, double newBeta);	
  
  inline double getBeta(int pos) { return beta[pos]; };
  
  inline vector<double> getBetaVec() { return beta; };
  
  inline int getBetaSize() { return beta.size(); };
  
  inline double getHessian(int pos) { return diagXTWX[pos]; }
  
  inline int getn() const { return X.n; }
  
  inline int getp() const { return X.p; }
  
  inline void activate(int pos) {};
  
  vector<double> getDerivativeVec();
  
  bool isExtreme() { return false; };
};

// the class that keeps the data for calculating the derivative
// intended for logistic regression
// only differs from the standard class 
class QuadraticDerivativeLogistic : public QuadraticDerivative {
public:
  int n; // number of rows in X
  int p; // number of cols in X
  
  SparseMatrix X; // copy of vector X
  vector<double> y; // copy of vector y
  vector<double> z; // adjusted response
  vector<double> w; // observation weights
  vector<double> beta; // copy of R object that stores the current value of beta
  SparseMatrix WX; // stores a sparse matrix of WX
  vector<double> Xbeta; // stores the current value of Xbeta
  vector<double> diagXTWX;
  vector<double> XTWz;
  
  double cutoff;
  bool probIsNaN;
public:
  /*
   Returns the value of the quadratic derivative
   Encapsulates all necessary pre-computation in the Constructor
   Assumes that the dimensions of the objects are correct
   */
  QuadraticDerivativeLogistic() {};
  
  // Constructor
  QuadraticDerivativeLogistic(const SparseMatrix& X, const vector<double>& y, const vector<double>& w, const vector<double>& beta);

  // make the appropriate calculations in the constructor after the data has been set
  void constructorWorker();
  
  ~QuadraticDerivativeLogistic() {};
  
  // calculate the derivative for the position pos
  double getDerivative(int pos);	
  
  inline double getProb(const int pos) const;
  
  void updateBeta(int pos, double newBeta);	
  
  inline double getBeta(int pos) { return beta[pos]; };
  
  inline vector<double> getBetaVec() { return beta; };
  
  inline int getBetaSize() { return beta.size(); };
  
  inline double getHessian(int pos) { return diagXTWX[pos]; }
  
  inline int getn() const { return X.n; }
  
  inline int getp() const { return X.p; }
  
  inline void activate(int pos) {};
  
  vector<double> getDerivativeVec();
  
  bool isExtreme();
  
};

