#include "FusedLasso_utils.h"

using namespace std;

double RelDif(double a, double b)
{
  double c = std::abs(a);
  double d = std::abs(b);
  d = std::max(c, d);
  return d == 0.0 ? 0.0 : std::abs(a - b) / d;
}

double RelDifNoAbs(double a, double b)
{
  double c = std::abs(a);
  double d = std::abs(b);
  d = std::max(c, d);
  return d == 0.0 ? 0.0 : (a - b) / d;
}

double maxDiffDoubleVec(const vector<double>& x, const vector<double>& y) {
  double maxDiff = 0;
  if(x.size() != y.size()) {
    return infinite;
  }
  for(size_t i = 0; i < x.size(); ++i) {
    double diff = std::abs(x[i] - y[i]) ;
    if(diff > maxDiff) maxDiff = diff ;
  }
  return maxDiff;
}


int numNonZero(const vector<double>& x) {
  int nonZero = 0;
  for(const auto& xi : x) if(xi != 0) nonZero++ ;
  return nonZero;
}


void printVector(const vector<double>& x, ostream& outStream) {
  for(size_t i = 0; i < x.size(); ++i) outStream << i << ":" << x[i] << ", ";
  outStream << endl;
}

void printVector(const vector<int>& x, ostream& outStream) {
  for(size_t i = 0; i < x.size(); ++i) outStream << i << ":" << x[i] << ", ";
  outStream << endl;
}

void printList(const list<double>& x, ostream& outStream) {
  list<double>::const_iterator listIt; 
  for(listIt = x.begin(); listIt != x.end(); ++listIt) {
    outStream << *listIt << ", ";
  }
  outStream << endl;
}

void printList(const list<int>& x, ostream& outStream) {
  list<int>::const_iterator listIt; 
  for(listIt = x.begin(); listIt != x.end(); ++listIt) {
    outStream << *listIt << ", ";
  }
  outStream << endl;
}


void printMatrix(vector<double>& X, int n, int p, ostream& outStream) {
  outStream << "n: " << n << "   p: " << p << endl;
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < p; ++j) {
      outStream << X[j*n + i] << ", ";
    }
    outStream << endl;
  }
}

namespace Steps {

void findRootL2(vector<double>& beta, double deriv, const double hessian, const int pos, const double zeroStepSize, const vector<int>& l2Idx, const vector<double>& l2Mult) {
  
  // first gives the position of the hessianchange, second the amoung of hessianChange
  
  double curDeriv = deriv - beta[pos] * hessian;
  double curHessian = hessian;
  for(size_t i = 0; i < l2Idx.size(); ++i) {
    curDeriv += l2Mult[i] * 2 * (0 - beta[l2Idx[i]]);
    curHessian += l2Mult[i] * 2;
  }
  
  if(curDeriv > zeroStepSize) {
    beta[pos] = - (curDeriv-zeroStepSize)/ curHessian;
  }
  else if(curDeriv < -zeroStepSize) {
    beta[pos] = - (curDeriv + zeroStepSize) / curHessian;
  }
  else {
    beta[pos] = 0;
  }
  
  return;
}


// finding the optimum using the huber penalty; huberParam is the scaling factor for the
// Huber penalty (the larger the more similar to L1 penalty)
void findRootHuber(vector<double>& beta, double deriv, const double hessian, const int pos, const double zeroStepSize, const vector<int>& huberIdx, const vector<double>& huberMult, double huberParam, vector<pair<double,double>>& scratch) {
  double doubleMax = numeric_limits<double>::max();
  scratch.resize(2*huberIdx.size()+1);
  vector<pair<double, double>>& huberPosVec = scratch;
  double derivSum = zeroStepSize;
  int counter = 0;
  for(size_t i = 0; i < huberIdx.size(); ++i) {
    huberPosVec[counter] = pair<double, double>(beta[huberIdx[i]] - 1/huberParam, huberMult[i] * huberParam);
    counter++;
    huberPosVec[counter] = pair<double, double>(beta[huberIdx[i]] + 1/huberParam, - huberMult[i] * huberParam);
    counter++;
    derivSum += huberMult[i];
  }
  huberPosVec[counter] = pair<double, double>(0, doubleMax);
  
  // sort them
  sort(huberPosVec.begin(), huberPosVec.begin() + huberPosVec.size());
  
  double curBeta = huberPosVec[0].first - 1 / huberParam;
  double curDeriv = deriv + (curBeta - beta[pos]) * hessian - derivSum;
  //    cout << "deriv: " << deriv << "DerivSum " << derivSum << "hessian " << "beta[pos]" << beta[pos] << "curBeta" << curBeta << "curderiv" << curDeriv <<endl;;
  
  
  double curHessian = hessian;
  
  //    cout << "CurBeta: " << curBeta << " curDeriv: " << curDeriv << " curHessian: " << curHessian << endl;
  
  if(curDeriv > 0) { // the root is to the left
    beta[pos] = curBeta - curDeriv/curHessian;
    return;
  }
  
  for(size_t i = 0; i < huberPosVec.size(); ++i) {
    //        cout << "pos: " << huberPosVec[i].first << " val: " << huberPosVec[i].second << endl;
    curDeriv += curHessian * (huberPosVec[i].first - curBeta);
    curBeta = huberPosVec[i].first;
    //        cout << "beta: " << curBeta << "deriv " << curDeriv;
    if(curDeriv > 0) { // found the root
      beta[pos] = curBeta - curDeriv/curHessian;
      //    cout << "CurBeta: " << curBeta << " curDeriv: " << curDeriv << " curHessian: " << curHessian << endl;
      return;
    }
    // check if it is special
    if(huberPosVec[i].second == doubleMax) {
      
      curDeriv += 2 * zeroStepSize;
      //    cout << "CurBeta: " << curBeta << " curDeriv: " << curDeriv << " curHessian: " << curHessian << endl;
      if(curDeriv > 0) {
        beta[pos] = huberPosVec[i].first;
        return;
      }
    }
    else {
      // regular huber penalty
      curHessian += huberPosVec[i].second;
      //    cout << "CurBeta: " << curBeta << " curDeriv: " << curDeriv << " curHessian: " << curHessian << endl;
    }
  }
  
  // still not found 
  beta[pos] = curBeta - curDeriv/curHessian;
  return;
}

void findRootL1(vector<double>& beta, double deriv, const double hessian, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepMult, vector<pair<double,double>>& scratch) {
  scratch.resize(stepIdx.size()+1);
  vector<pair<double, double>>& stepPosVec = scratch;
  double derivSum = zeroStepSize;
  for(size_t i = 0; i < stepIdx.size(); ++i) {
    stepPosVec[i] = pair<double, double>(beta[stepIdx[i]], 2 * stepMult[i]);
    derivSum += stepMult[i];
  }
  stepPosVec[stepIdx.size()] = pair<double, double>(0, 2 * zeroStepSize);
  
  // sort them
  sort(stepPosVec.begin(), stepPosVec.begin() + stepPosVec.size());
  
  double curBeta = stepPosVec[0].first;
  double curDeriv = deriv + (curBeta - beta[pos]) * hessian - derivSum;
  
  
  if(curDeriv > 0) { // the root is to the left
    if(hessian > 1e-8) {
      beta[pos] = curBeta - curDeriv/hessian;
    }
    else {
      beta[pos] = curBeta; // do nothing as there is a problem but we do not want to crash
    }
    return;
  }
  
  for(size_t i = 0; i < stepPosVec.size(); ++i) {
    curDeriv += hessian * (stepPosVec[i].first - curBeta);
    curBeta = stepPosVec[i].first;
    
    // regular huber penalty
    if(curDeriv > 0) { // found the root; must have positive hessian
      beta[pos] = curBeta - curDeriv/hessian;
      return;
    }
    // now add the step
    curDeriv += stepPosVec[i].second;
    if(curDeriv > 0) {
      beta[pos] = curBeta;
      return;
    }
  }
  
  // still not found 
  if(hessian > 1e-8) {
    beta[pos] = curBeta - curDeriv/hessian;
  }
  else {
    beta[pos] = curBeta;
  }
  return;
}


/*
  * Update beta, return true when finished (moves at most over one step)
*/
  static bool updateBeta(vector<double>& beta, const int pos, double& deriv, const double& slope, const double& zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize) {
    double curDeriv = deriv;
    double curDerivDelta = 0;
    double nextSmallerPos = -numeric_limits<double>::max();
    double nextSmallerSize = 0;
    double nextLargerPos = numeric_limits<double>::max();
    double nextLargerSize = 0;
    double curBeta = beta[pos];
    
    if(curBeta > 0) {
      curDeriv += zeroStepSize;
      nextSmallerPos = 0;
      nextSmallerSize = zeroStepSize;
    }
    else if(curBeta < 0) {
      curDeriv -= zeroStepSize;
      nextLargerPos = 0;
      nextLargerSize = zeroStepSize;
    }
    else {
      curDerivDelta = zeroStepSize;
    }
    
    
    // calculate current derivative, also the next position that is smaller or larger
    for(size_t i = 0; i < stepIdx.size(); ++i) {
      if(beta[stepIdx[i]] > curBeta) {
        curDeriv -= stepSize[i];
        if(beta[stepIdx[i]] < nextLargerPos) {
          nextLargerPos = beta[stepIdx[i]];
          nextLargerSize = stepSize[i];
        }
      }
      else if(beta[stepIdx[i]] < curBeta) {
        curDeriv += stepSize[i];
        if(beta[stepIdx[i]] > nextSmallerPos) {
          nextSmallerPos = beta[stepIdx[i]];
          nextSmallerSize = stepSize[i];
        }
      }            
      else { // both are equal
        curDerivDelta += stepSize[i];
      }
    }
    
    // check if we are on a step and optimal
    if(curDeriv <= curDerivDelta && curDeriv >= -curDerivDelta) {
      return(true);
    }
    
    // cout << "CurBeta: " << curBeta << endl;
    // cout << "CurDeriv: " << curDeriv << endl;
    // cout << "nextLarger: " << nextLargerPos << endl;
    // cout << "nextSmaller: " << nextSmallerPos << endl;
    
    // now check if a step to the left or right should be made
    double newBeta;
    double newStepSize = tolerance;
    if(curDeriv > 0) {
      curDeriv -= curDerivDelta; // adjust for being on a step
      newBeta = curBeta - curDeriv / slope;
      if(newBeta <= nextSmallerPos) { // at least at the left step
        newBeta = nextSmallerPos;
        newStepSize = nextSmallerSize;
      }
    }
    else {
      curDeriv += curDerivDelta;
      newBeta = curBeta - curDeriv / slope;
      if(newBeta >= nextLargerPos) { // at least at the right step
        newBeta = nextLargerPos;
        newStepSize = nextLargerSize;
      }
    }
    
    // check if it is below the the smaller or above the larger
    // cout << "NewBeta: " << newBeta << endl;
    // cout << "curBeta: " << curBeta << endl;
    // cout << "Slope: " << slope << endl;
    
    deriv += (newBeta - curBeta) * slope;
    beta[pos] = newBeta;
    if(curDeriv <= newStepSize && curDeriv >= -newStepSize) {
      return(true);
    }
    else
    {
      return(false);
    }
  }


/*
  * Includes the step given into the calculation and updates steps on the left or right if necessary
*/
  
  void findRoot(vector<double> &beta, double derivQuad, const double slope, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize) {
    
    // catch if the slope is 0; this should only occur if the derivative is also 0;
    // then set the beta to 0
    if(fabs(slope) < 1e-10) {
      if(fabs(derivQuad) < 1e-10) {
        beta[pos] = 0;
        return;
      }
      else {
        Rcpp::stop("There has been a problem with a variable in Steps with 0 slope but non-zero derivative\n derivQuad %f slope: %f\n", derivQuad, slope, "\n");
      }
    }
    
    // generate a list with the steps
    while(!updateBeta(beta, pos, derivQuad, slope, zeroStepSize, stepIdx, stepSize)) {};
  }

} // namespace Steps

