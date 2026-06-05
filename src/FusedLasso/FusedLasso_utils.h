#pragma once

#include "RcppArmadillo.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <limits>
#include <list>
#include <map>
#include <queue>
#include <cmath>

using std::vector;
using std::list;
using std::map;
using std::queue;
using std::ostream;
using std::pair;

const double tolerance =1.0e-8;
const double infinite = std::numeric_limits<double>::max();
const int infiniteInt = std::numeric_limits<int>::max();

double RelDif(double a, double b);
double RelDifNoAbs(double a, double b);

double maxDiffDoubleVec(const vector<double>&x, const vector<double>& y);

inline int signum(double x) {return((x>0)-(x<0));};

int numNonZero(const vector<double>& x);

void printVector(const vector<double>& x, ostream& outStream);

void printVector(const vector<int>& x, ostream& outStream); 

void printList(const list<double>& x, ostream& outStream); 

void printList(const list<int>& x, ostream& outStream); 

void printMatrix(vector<double>& X, int n, int p, ostream& outStream);

/* NAMESPACE STEPS
 * Given the current betas, it finds the root of the function
 */

namespace Steps {
  void findRoot(vector<double> &beta, double derivQuad, const double slope, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
  void findRootL1(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize, vector<pair<double,double>>& scratch);
  void findRootL2(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& l2Idx, const vector<double>& l2Mult);
  void findRootHuber(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& huberIdx, const vector<double>& huberMult, double huberParam, vector<pair<double,double>>& scratch);
}

// C++11 smart pointers
#include <memory>
using std::shared_ptr;
using std::make_shared;

