#pragma once

#include "FusedLasso_utils.h"
#include "FusedLasso_enums.h"
#include "FusedLasso_optim.h"
#include "FusedLasso_data_struct.h"


// main class of the program; stores important data and administrate other important classes
class FusedLasso {
  
  // FIELDS
  
public:
  SparseMatrix X;
  vector<double> y;
  vector<double> wObs;
  vector<double> beta;
  
  vector<int> fusions; // stores how the beta are currently fused
  vector<int> fusedGroupSize; // how large is each of the groups
  
  vector<int> lastDetailedFusions; //stores the result of the last fusion
  vector<double> wLambda1;

  // Cache for the fused matrix and graph connections — invalidated when fusions change
  SparseMatrix cachedFusedX_;
  vector<vector<int>>  cachedConnectionsFused_;
  vector<vector<double>> cachedWLambda2Fused_;
  vector<double> cachedWLambda1Coordinate_;
  bool fusedCacheValid_ = false;

  void rebuildFusedCache_();
  
  // quadratic derivatives under the current fusions
  shared_ptr<QuadraticDerivative> quadratic;
  
  Graph pg;
  
  int n;
  int p;
  double lambda1;
  double lambda2;
  
  // these are used to initialize flc
  double accuracy;
  int maxIterInner;
  int maxIterOuter;
  int innerIterNum;
  int outerIterNum;
  int maxActivateVars;
  int fusionLevel;
  
  regEnum regType;
  
  // METHODS
  
  // calculate the pull on every variable (not taking into account 
  // variables that are on the same level i.e. could be fused)
  vector<double> getPulls();
  
  // find the fusions from the univariate tensions
  void getEqualFusions(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool zeroSingle=true, double accFactor=1);
  void getSplitFusionsActive(vector<int>& newFusions, vector<int>& newFusedGroupSize);
  void getSplitFusionsInactive(vector<int>& newFusions, vector<int>& newFusedGroupSize);

  // shared implementation for getSplitFusionsActive / getSplitFusionsInactive
  void buildFusionGroups_(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool splitZeroNodes);
  
  // check if two fusions are equal 
  bool areFusionsEqual(vector<int> &fusion1, vector<int> &fusion2);

  // adjusts the nodePull so that the maximum-flow graph can be run on it
  void makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2);
  
  // given nodes and groups, adds the complementary groups as defined by the given graph
  void addComplementaryGroups(vector<vector<int>>& groups, vector<int>& nodes);
  // compute the mean of the current group
  double calcGroupAverage(int pos, vector<double>& nodePull, vector<double>& beta);

  void sortAllGroups(vector<vector<int>>& x);

  // Constructor
  // Here regType can be GAUSSIAN or LOGISTIC

  // construct with sparse matrix as input
  FusedLasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType=regEnum::GAUSSIAN);
  
  ~FusedLasso();
  
  // translate the beta to a fused beta using the fusions saved in the class
  vector<double> translateOriginalToFused();
  
  // find the new fusions and set them
  // returns false if they are the same as the old ones
  // i.e. the algorithm has converged
  bool identifyNewFusions(double lastIterChange, FusionStrategy maxFusion);
  bool identifyNewFusionsHuber();
  
  // initialize the FusedLassoCoordinate object
  // run it to convergence and update all the objects in the current class
  // to the new values of beta
  bool runFused(penEnum penType);
  
  // run the fused lasso for survival data
  bool runFusedGeneral(penEnum penType);
  
  // After initialization, run the whole algorithm
  bool runAlgorithm(FusionStrategy maxFusion);
  bool runAlgorithmHuber();
  bool runAlgorithmL2();

  SparseMatrix runAlgorithm(penEnum penType, FusionStrategy maxFusion, const vector<double>& lambda1Vec, const vector<double>& lambda2Vec, const int maxNonZero, vector<bool>& success, vector<int>& outerIterNumVec, vector<int>& innerIterNumVec, bool verbose=false);
  
  // set new values for lambda1 and lambda2
  void setNewLambdas(double lambda1, double lambda2);
  
  inline int getOuterIterNum() {return outerIterNum ;} ;
  inline int getInnerIterNum() {return innerIterNum ;};
  
  // finds the maximum lambda1 at which the first non-exempt variable enters the model
  double findMaxLambda1(const vector<int>& exemptVars);
  
  // finds the value of lambda2 at which the first group of variables breaks apart
  double findMaxLambda2(double lambda1);
  
  vector<double> getBeta() { return beta; };
  
  // adds a column to betaMat with the new solution
  void getBeta(SparseMatrix& betaMat);
  
  void printBetaAndGroup(ostream& outStream);
  void printGroups(vector<vector<int>> x, ostream& outStream);
  void checkSolution(ostream& outStream);
};


