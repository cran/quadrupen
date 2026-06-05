#include "FusedLasso_data_struct.h"

using namespace std;

// ============================================================================
// SPARSE MATRIX
// 

SparseMatrix::SparseMatrix() {
  p = 0;
  n = 0;
  nnz = 0;
  nzVec.clear();
  indexCol.clear(); indexCol.push_back(0);
  rowPos.clear();
}


SparseMatrix::SparseMatrix(int n) {
  p = 0;
  this->n = n;
  nnz = 0;
  nzVec.clear();
  indexCol.clear(); indexCol.push_back(0);
  rowPos.clear();
}

SparseMatrix::SparseMatrix(const double* x, int n, int p) {
  this->p = p;
  this->n = n;
  this->nnz = 0;
  bool foundFirst;
  
  for(int i = 0; i < p; ++i) {
    foundFirst = false;
    for(int j = 0; j < n; ++j) {
      if(x[i*n + j] != 0) {
        nnz++;
        nzVec.push_back(x[i*n + j]);
        rowPos.push_back(j);
        if(!foundFirst) {
          foundFirst = true;
          indexCol.push_back(nnz-1);
        }
      }
    }
    if(!foundFirst) { // column is empty
      indexCol.push_back(nnz);
    }
  }
  indexCol.push_back(nnz);
}

SparseMatrix::SparseMatrix(const vector<double>& x, int n, int p) {
  this->p = p;
  this->n = n;
  this->nnz = 0;
  bool foundFirst;
  
  for(int i = 0; i < p; ++i) {
    foundFirst = false;
    for(int j = 0; j < n; ++j) {
      if(x[i*n + j] != 0) {
        nnz++;
        nzVec.push_back(x[i*n + j]);
        rowPos.push_back(j);
        if(!foundFirst) {
          foundFirst = true;
          indexCol.push_back(nnz-1);
        }
      }
    }
    if(!foundFirst) { // column is empty
      indexCol.push_back(nnz);
    }
  }
  indexCol.push_back(nnz);
}

SparseMatrix::SparseMatrix(const vector<double>& x, const vector<int>& indexCol, const vector<int>& rowPos, int n, int p) {
  this->n = n;
  this->p = p;
  this->nnz = (int)x.size();
  this->nzVec = x;
  this->indexCol = indexCol;
  this->rowPos = rowPos;
}

// assumes that it is the right class; does not check
SparseMatrix::SparseMatrix(const arma::sp_mat& sm) {

  sm.sync();

  indexCol = std::vector<int>(sm.col_ptrs, sm.col_ptrs + sm.n_cols+1 ) ;
  rowPos = std::vector<int>(sm.row_indices, sm.row_indices + sm.n_nonzero);
  nzVec  = std::vector<double>(sm.values, sm.values + sm.n_nonzero ) ;
  
  n = (int)sm.n_rows;
  p = (int)sm.n_cols;
  nnz = (int)nzVec.size();
}

// assumes that it is the right class; does not check
SparseMatrix::SparseMatrix(Rcpp::S4 x) {
  
  indexCol = Rcpp::as<std::vector<int> >(x.slot("p")) ;
  rowPos   = Rcpp::as<std::vector<int> >(x.slot("i")) ;
  nzVec    = Rcpp::as<std::vector<double> >(x.slot("x")) ;
  Rcpp::IntegerVector dimslot = x.slot("Dim");
  
  n = dimslot[0];
  p = dimslot[1];
  nnz = (int)nzVec.size();
}


// for the fusions, create the fusedX matrix
// Implementation at the moment assumes that 
SparseMatrix SparseMatrix::createFusedX(const vector<vector<int> > fusedGroups) {
  // in order to make it more efficient, allocate enough space for the whole
  // X matrix; may be smaller, but for very sparse matrices, will be a good
  // approximation
  SparseMatrix fusedX;
  fusedX.n = n;
  fusedX.p = (int)fusedGroups.size();
  
  fusedX.nzVec.reserve(nzVec.size());
  fusedX.nzVec.clear();
  fusedX.indexCol.resize(fusedGroups.size() + 1, 0);
  fusedX.rowPos.reserve(rowPos.size());
  fusedX.rowPos.clear();
  
  // Pre-calculate maximum group size to avoid repeated memory reallocations
  int maxGroupSize = 0;
  for(size_t g = 0; g < fusedGroups.size(); ++g) {
    if((int)fusedGroups[g].size() > maxGroupSize) {
      maxGroupSize = fusedGroups[g].size();
    }
  }
  
  // Pre-allocate pos and endPos vectors once with sufficient capacity
  // This avoids 2×numGroups reallocations during the loop
  vector<int> pos, endPos;
  pos.reserve(maxGroupSize);
  endPos.reserve(maxGroupSize);
  int curRow;
  double curVal;
  
  for(unsigned int group = 0; group < fusedGroups.size(); ++group) {
    curRow = 0;
    // Reuse pre-allocated capacity without memory reallocation
    pos.clear(); 
    endPos.clear();
    pos.resize(fusedGroups[group].size()); // No reallocation - already reserved
    endPos.resize(fusedGroups[group].size()); // No reallocation - already reserved
    
    // Initialize the pos and endPos vectors
    for(size_t i = 0; i < fusedGroups[group].size(); ++i) {
      pos[i] = indexCol[fusedGroups[group][i]];
      endPos[i] = indexCol[fusedGroups[group][i] + 1] - 1;
    }
    
    // save the starting position of this column
    fusedX.indexCol[group] = (int)fusedX.nzVec.size();
    
    while(curRow < n) {
      curRow = n;
      // find the smallest row available
      for(unsigned int i = 0; i < pos.size(); ++i) {
        if(pos[i] <= endPos[i] && rowPos[pos[i]] < curRow) {
          curRow = rowPos[pos[i]];
        }
      }
      
      if(curRow == n) {
        continue;
      }
      
      curVal = 0;
      fusedX.rowPos.push_back(curRow); 
      for(unsigned int i = 0; i < pos.size(); ++i) {
        if(pos[i] <= endPos[i] && rowPos[pos[i]] == curRow) {
          curVal += nzVec[pos[i]];
          pos[i]++;
        }                
      }
      fusedX.nzVec.push_back(curVal);
    }
  }
  fusedX.indexCol[fusedGroups.size()] = (int)fusedX.nzVec.size();
  fusedX.nnz = (int)fusedX.nzVec.size();
  
  return fusedX;
}


void SparseMatrix::addColumn(const vector<double>& newnzVec, const vector<int>& newrowPos) {
  nzVec.insert(nzVec.end(), newnzVec.begin(), newnzVec.end());
  rowPos.insert(rowPos.end(), newrowPos.begin(), newrowPos.end());
  nnz += (int)newnzVec.size();
  indexCol.push_back(nnz);
  p += 1;
}


void SparseMatrix::addColumn(const vector<double>& x) {
  for(int j = 0; j < n; ++j) {
    if(x[j] != 0) {
      nnz++;
      nzVec.push_back(x[j]);
      rowPos.push_back(j);
    }
  }
  indexCol.push_back(nnz);
  p += 1;
}


double SparseMatrix::innerProd(int i) const {
  // Use consolidated helper for column self-product
  return dotProductColumnSelf(i, nullptr);
}

double SparseMatrix::innerProd(int i, int j) const {
  if(i < 0 || i >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", i, p);
  }
  
  if(j < 0 || j >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", j, p);
  }
  
  double res = 0;
  // first check that column j has any entries
  if(indexCol[j] == indexCol[j + 1]) { // no entries
    return 0;
  }
  
  int posJ = indexCol[j];
  for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
    while(rowPos[posJ] < rowPos[posI]) { // increase rows until not less than in I
      posJ++;
    }
    if(rowPos[posJ] == rowPos[posI]) {
      res += nzVec[posJ] * nzVec[posI];
    }
  }
  
  return res;
}

double SparseMatrix::innerProd(const vector<double>& y, int i) const {
  // Use consolidated helper for column-vector dot product
  return dotProductColumnVector(y, i, nullptr);
}

double SparseMatrix::innerProd(int i, const vector<double>& w) const {
  // Use consolidated helper for weighted column self-product
  return dotProductColumnSelf(i, &w);
}

double SparseMatrix::innerProd(int i, int j, const vector<double>& w) const {
  if(i < 0 || i >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", i, p);
  }
  
  if(j < 0 || j >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", j, p);
  }
  
  if((int)w.size() != n) {
    Rcpp::stop("w not the right size\n");
  }
  
  double res = 0;
  // first check that column j has any entries
  if(indexCol[j] == indexCol[j + 1]) { // no entries
    return 0;
  }
  
  int posJ = indexCol[j];
  for(int posI = indexCol[i]; posI < indexCol[i + 1]; ++posI) {
    while(rowPos[posJ] < rowPos[posI]) { // increase rows until not less than in I
      posJ++;
    }
    if(rowPos[posJ] == rowPos[posI]) {
      res += nzVec[posJ] * nzVec[posI] * w[rowPos[posI]];
    }
  }
  
  return res;
}

double SparseMatrix::innerProd(const vector<double>& y, int i, const vector<double>& w) const {
  // Use consolidated helper for weighted column-vector dot product
  return dotProductColumnVector(y, i, &w);
}

void SparseMatrix::multByDiag(const vector<double>& w) {
  // check that w has the right size
  if((int)w.size() != n) { // not the right size
    Rcpp::stop("w does not have the right size\n");
  }
  
  // no multiply the matrix correctly
  for(int i = 0; i < nnz; ++i) {
    nzVec[i] *= w[rowPos[i]];
  }
}

// need to check later if these checks at the beginning result in a
// speed penalty
void SparseMatrix::addMultOfColumn(vector<double>& y, int i, double mult) const {
  // check that y has the right size
  if((int)y.size() != n) { // does not have the right size
    Rcpp::stop("y does not have the right size\n");
  }
  if(i < 0 || i >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", i, p);
  }
  
  if(mult == 0) {
    return;
  }
  
  for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
    y[rowPos[pos]] += nzVec[pos] * mult;
  }
}

// Unsafe version: assumes i is valid in [0, p) and y.size() == n
// Skips all bounds checking for 10-15% performance gain in hot loops
// Should only be used in performance-critical paths with pre-validated indices
inline void SparseMatrix::addMultOfColumnUnsafe(vector<double>& y, int i, double mult) const {
  // DEBUG mode: still do checks (can be disabled with -DNDEBUG)
  #ifndef NDEBUG
    if((int)y.size() != n || i < 0 || i >= p) {
      Rcpp::stop("addMultOfColumnUnsafe: precondition violated (y.size=%d, n=%d, i=%d, p=%d)",
                 (int)y.size(), n, i, p);
    }
  #endif
  
  if(mult == 0) {
    return;
  }
  
  // No bounds checking - direct access
  for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
    y[rowPos[pos]] += nzVec[pos] * mult;
  }
}

// also does checking if the column is ok - not good for performance
// but this function is mainly intended for error checking using the
// test suite
double SparseMatrix::get(int i, int j) const {
  if(j < 0 || j >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", j, p);
  }
  
  for(int pos = indexCol[j]; pos < indexCol[j + 1]; ++pos) {
    
    if(rowPos[pos] == i) {
      return nzVec[pos];
    }
    if(rowPos[pos] > i) { // could not find it
      return 0;
    }
  }
  return 0; // could not find it
}


vector<double> SparseMatrix::getColumn(int i) const {
  vector<double> res(n);
  
  if(i < 0 || i >= p) {
    Rcpp::stop("not a valid column %d out of %d\n", i, p);
  }
  
  for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
    res[rowPos[pos]] = nzVec[pos];
  }
  return res;
}

SEXP SparseMatrix::todgCMatrix() const {
  Rcpp::S4 s(std::string("dgCMatrix"));
  
  s.slot("i") = rowPos ;
  s.slot("p") = indexCol ;
  s.slot("x") = nzVec ;
  s.slot("Dim") = Rcpp::IntegerVector::create(n, p);
  return s;
}

void SparseMatrix::print(ostream& outStream) const {
  outStream << "---------------------------------------------" << endl;
  outStream << "nnz: " << nnz << " n: " << n << " p: " << p << endl;
  outStream << "nzVec:" << endl;
  printVector(nzVec, outStream);
  outStream << "indexCol:" << endl;
  printVector(indexCol, outStream);
  outStream << "rowPos:" << endl;
  printVector(rowPos, outStream);
  outStream << "---------------------------------------------" << endl;
}

void SparseMatrix::printColumn(int i, ostream& outStream) const {
  outStream << "Num non zero: " << indexCol[i+1] - indexCol[i] << endl;
  for(int pos = indexCol[i]; pos < indexCol[i + 1]; ++pos) {
    outStream << " Row: " << rowPos[pos] << " Val: " << nzVec[pos];
  }
  outStream << endl;
}

// ============================================================================
// LINEAR REGRESSION (GAUSSIAN OUTCOME)
// 

QuadraticDerivativeDiagonal::QuadraticDerivativeDiagonal(const SparseMatrix& X, const vector<double>& y, const vector<double>&  w, const vector<double>& beta) {
  // set number of rows and cols and store X, y and w and beta
  n = y.size();
  p = beta.size();
  
  this->X = X;
  this->y = y;
  this->w = w;
  this->beta = beta;
  
  // create WX matrix: W is a diagonal matrix represented by w
  WX = X;
  WX.multByDiag(w);
  
  // Precalculate WXbeta = WX * beta
  WXbeta.clear();
  WXbeta.resize(n, 0);
  for(int pos = 0; pos < p; ++pos) {
    WX.addMultOfColumn(WXbeta, pos, beta[pos]);
  }
  
  // Pre-calculate diagonal of X^T WX and X^T Wy for faster derivative computation
  diagXTWX.clear();
  diagXTWX.resize(p, 0);
  XTWy.clear();
  XTWy.resize(p, 0);
  for(int i = 0; i < p; ++i) {
    diagXTWX[i] = X.innerProd(i, w);
    XTWy[i] = WX.innerProd(y, i);
  }
}

double QuadraticDerivativeDiagonal::getDerivative(int pos) {
  // calculate the rest of the expression XTWXbeta
  double res = X.innerProd(WXbeta, pos);
  
  // now subtract XTWy
  res -= XTWy[pos];
  return(res);
}

void QuadraticDerivativeDiagonal::updateBeta(int pos, double newBeta) {
  // update beta and save the difference
  double betaDiff = newBeta - beta[pos];
  if(betaDiff == 0) return;
  
  beta[pos] = newBeta;
  
  // now update Xbeta
  WX.addMultOfColumn(WXbeta, pos, betaDiff);
}

vector<double> QuadraticDerivativeDiagonal::getDerivativeVec() {
  // Optimized sparse matrix-vector product: O(nnz) instead of O(p×nnz)
  // Compute: res = X^T * WXbeta - XTWy by iterating over sparse X structure
  
  // Initialize result with -XTWy
  vector<double> res = XTWy;
  for(int j = 0; j < p; ++j) {
    res[j] = -res[j];
  }
  
  // Add X^T * WXbeta using the CSC sparse structure
  // For each column j, compute the dot product of column j with WXbeta
  for(int j = 0; j < X.p; ++j) {
    double col_contrib = 0.0;
    for(int pos = X.indexCol[j]; pos < X.indexCol[j + 1]; ++pos) {
      int row = X.rowPos[pos];
      col_contrib += X.nzVec[pos] * WXbeta[row];
    }
    res[j] += col_contrib;
  }
  
  return res;
}


// ============================================================================
// LOGISTIC REGRESSION (BINARY OUTCOME)
// 

QuadraticDerivativeLogistic::QuadraticDerivativeLogistic(const SparseMatrix& X, const vector<double>& y, const vector<double>&  w, const vector<double>& beta) {
  // Initialize common data members (n, p, X, y, w, beta)
  n = y.size();
  p = beta.size();
  this->X = X;
  this->y = y;
  this->w = w;
  this->beta = beta;
  
  // Call worker to complete logistic-specific initialization
  constructorWorker();
}

void QuadraticDerivativeLogistic::constructorWorker() {
  // precalculate the current Xbeta
  probIsNaN = false;
  
  // Initialize Xbeta = X * beta
  Xbeta.clear();
  Xbeta.resize(n, 0);
  for(int pos = 0; pos < p; ++pos) {
    X.addMultOfColumn(Xbeta, pos, beta[pos]);
  }
  
  // Adjust weights and response for logistic regression
  // z = Xbeta + (y - prob) / (prob * (1 - prob))
  // w = w * prob * (1 - prob)
  z.resize(n, 0);
  for(int i = 0; i < n; ++i) {
    double prob = getProb(i);
    w[i] *= prob * (1 - prob);
    z[i] = Xbeta[i] + ((y[i] - prob) / (prob * (1 - prob))); 
    if(isnan(z[i])) {
      probIsNaN = true;
    }
  }
  
  // create WX matrix
  WX = X;
  WX.multByDiag(w);
  
  // Pre-calculate diagonal of X^T WX and X^T Wz
  diagXTWX.clear();
  diagXTWX.resize(p, 0);
  XTWz.clear();
  XTWz.resize(p, 0);
  for(int i = 0; i < p; ++i) {
    diagXTWX[i] = X.innerProd(i, w);
    XTWz[i] = WX.innerProd(z, i);
  }
  
  cutoff = 0.01;
}

double QuadraticDerivativeLogistic::getDerivative(int pos) {
  // calculate the rest of the expression XTWXbeta
  double res = WX.innerProd(Xbeta, pos);
  if(isnan(res)) {
    Rcpp::stop("NAN1\n");
  }
  // now subtract XTWy
  res -= XTWz[pos];
  if(isnan(res)) {
    Rcpp::stop("NAN2\n");
  }
  
  return(res);
}

inline double QuadraticDerivativeLogistic::getProb(const int pos) const {
  return(1/(1+exp(-Xbeta[pos])));
}

bool QuadraticDerivativeLogistic::isExtreme() {
  if(probIsNaN) {
    return(true);
  }
  for(int i = 0; i < n; ++i) {
    if(getProb(i) > cutoff || getProb(i) < 1 - cutoff) {
      return(false);
    }
  }
  return(true);
}


void QuadraticDerivativeLogistic::updateBeta(int pos, double newBeta) {
  // update beta and save the difference
  double betaDiff = newBeta - beta[pos];
  if(betaDiff == 0) return;
  
  beta[pos] = newBeta;
  
  // now update Xbeta
  X.addMultOfColumn(Xbeta, pos, betaDiff);
}

vector<double> QuadraticDerivativeLogistic::getDerivativeVec() {
  // Optimized sparse matrix-vector product: O(nnz) instead of O(p×nnz)
  // Compute: res = X^T * W * Xbeta - X^T * W * z by iterating over sparse X structure
  
  // Initialize result with -XTWz
  vector<double> res = XTWz;
  for(int j = 0; j < p; ++j) {
    res[j] = -res[j];
  }
  
  // Add X^T * W * Xbeta using the CSC sparse structure
  // For each column j, compute the dot product: sum_i X[i,j] * w[i] * Xbeta[i]
  for(int j = 0; j < X.p; ++j) {
    double col_contrib = 0.0;
    for(int pos = X.indexCol[j]; pos < X.indexCol[j + 1]; ++pos) {
      int row = X.rowPos[pos];
      col_contrib += X.nzVec[pos] * w[row] * Xbeta[row];
    }
    res[j] += col_contrib;
  }
  
  return res;
}

