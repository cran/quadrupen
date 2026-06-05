#include "FusedLasso_class.h"

using namespace std;

FusedLasso::FusedLasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType) {

  this->X = X;
  this->y = y;
  this->wObs = wObs;
  this->beta = beta;
  this->n = X.n;
  this->p = X.p;
  this->wLambda1 = wLambda1;
  this->regType = regType;
  this->pg = graph;
  this->lambda1 = lambda1;
  this->lambda2 = lambda2;

  this->maxIterInner = maxIterInner;
  this->maxIterOuter = maxIterOuter;
  this->accuracy = accuracy;
  this->maxActivateVars = maxActivateVars;

  this->fusedGroupSize.resize(p, 1);
  this->fusions.resize(p);
  std::iota(this->fusions.begin(), this->fusions.end(), 0);

  switch(regType) {
  case regEnum::GAUSSIAN:
    quadratic = make_shared<QuadraticDerivativeDiagonal>(this->X, this->y, this->wObs, this->beta);
    break;
  case regEnum::BINOMIAL:
    quadratic = make_shared<QuadraticDerivativeLogistic>(this->X, this->y, this->wObs, this->beta);
  }

  lastDetailedFusions = fusions;
  fusionLevel = 0;
}

FusedLasso::~FusedLasso() {}

vector<double> FusedLasso::getPulls() {
  vector<double> pulls(beta.size());
  for(size_t i = 0; i < pulls.size(); ++i) pulls[i] = quadratic->getDerivative(i);
  return pulls;
}

void FusedLasso::getEqualFusions(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool zeroSingle, double accFactor) {
  newFusions.clear(); newFusedGroupSize.clear();
  newFusions.resize(beta.size(), -1);
  int groupNum = 0;
  double accHere = accuracy * accFactor;

  for(int i = 0; i < (int)beta.size(); ++i) {
    if(newFusions[i] == -1) {
      vector<int> nodes = pg.connectedWithSameValue(i, beta, accHere);
      if(zeroSingle && beta[i] == 0) {
        for(int node : nodes) {
          newFusions[node] = groupNum++;
          newFusedGroupSize.push_back(1);
        }
      }
      else {
        newFusedGroupSize.push_back(nodes.size());
        for(int node : nodes) newFusions[node] = groupNum;
        ++groupNum;
      }
    }
  }
}

bool FusedLasso::areFusionsEqual(vector<int>& fusion1, vector<int>& fusion2) {
  if(fusion1.size() != fusion2.size()) return false;
  for(size_t i = 0; i < fusion1.size(); ++i) {
    if(fusion1[i] != fusion2[i]) return(false);
  }
  return true;
}

// Shared implementation for getSplitFusionsActive and getSplitFusionsInactive.
// splitZeroNodes=true: zero-beta nodes become individual singletons (active variant).
// splitZeroNodes=false: zero-beta nodes are split via max-flow (inactive variant).
void FusedLasso::buildFusionGroups_(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool splitZeroNodes) {
  newFusions.clear(); newFusedGroupSize.clear();
  newFusions.resize(p, -1);

  vector<double> nodePull = getPulls();
  makePullAdjustment(beta, nodePull, lambda2);
  vector<vector<int>> groups;
  int newNumFused = 0;

  for(int i = 0; i < p; ++i) {
    if(newFusions[i] == -1) {
      vector<int> nodes = pg.connectedWithSameValue(i, beta, accuracy);
      if(beta[i] == 0) {
        if(splitZeroNodes) {
          for(int node : nodes) {
            beta[node] = 0;
            newFusedGroupSize.push_back(1);
            newFusions[node] = newNumFused++;
          }
          continue;
        } else {
          groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
          vector<vector<int>> positiveGroups = pg.splitGroup(nodes, nodePull, lambda1, lambda2, true);
          groups.insert(groups.end(), positiveGroups.begin(), positiveGroups.end());
        }
      } else if(splitZeroNodes) {
        // Active variant: split non-zero groups via max-flow
        if(beta[i] > 0) groups = pg.splitGroup(nodes, nodePull,  lambda1, lambda2, false);
        else             groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
      } else {
        // Inactive variant: keep non-zero groups intact
        groups.assign(1, nodes);
      }

      try { addComplementaryGroups(groups, nodes); }
      catch(std::exception& ex) { forward_exception_to_r(ex); }

      for(size_t j = 0; j < groups.size(); ++j) {
        newFusedGroupSize.push_back(groups[j].size());
        for(int node : groups[j]) newFusions[node] = newNumFused;
        newNumFused++;
      }
    }
  }
}

void FusedLasso::getSplitFusionsActive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {
  buildFusionGroups_(newFusions, newFusedGroupSize, /*splitZeroNodes=*/true);
}

void FusedLasso::getSplitFusionsInactive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {
  buildFusionGroups_(newFusions, newFusedGroupSize, /*splitZeroNodes=*/false);
}


void FusedLasso::sortAllGroups(vector<vector<int>>& x) {
  for(auto& group : x) sort(group.begin(), group.end());
}

static bool groupComp(const vector<int>& a, const vector<int>& b) {
  return a.front() < b.front();
}

void FusedLasso::addComplementaryGroups(vector<vector<int>>& groups, vector<int>& nodes) {
  // Merge all nodes from groups, check for duplicates, then find complement in nodes
  vector<int> groupsMerged;
  for(auto& g : groups) groupsMerged.insert(groupsMerged.end(), g.begin(), g.end());

  sort(groupsMerged.begin(), groupsMerged.end());
  size_t before = groupsMerged.size();
  groupsMerged.erase(unique(groupsMerged.begin(), groupsMerged.end()), groupsMerged.end());
  if(before > groupsMerged.size()) throw std::out_of_range("A node was grouped twice...");

  vector<int> complNodes = pg.getComplement(groupsMerged, nodes);

  vector<vector<int>> complementGroups = pg.identifyConnectedGroups(complNodes);
  groups.insert(groups.end(), complementGroups.begin(), complementGroups.end());
  sortAllGroups(groups);
  sort(groups.begin(), groups.end(), groupComp);
}

vector<double> FusedLasso::translateOriginalToFused() {
  vector<double> fusedBeta(fusedGroupSize.size());
  for(size_t i = 0; i < fusions.size(); ++i) fusedBeta[fusions[i]] = beta[i];
  return fusedBeta;
}

bool FusedLasso::identifyNewFusionsHuber() {
  vector<int> newFusions;
  vector<int> newFusedGroupSize;

  // huberParam chosen large enough to approximate L1 closely during the pre-smoothing step
  const double huberParam = 1000;

  vector<int> singleFusions(p);
  std::iota(singleFusions.begin(), singleFusions.end(), 0);

  vector<vector<int>> connectionsSingle;
  vector<vector<double>> wLambda2Single;
  pg.getFusedConnectionsWeights(singleFusions, p, connectionsSingle, wLambda2Single);

  FusedLassoCoordinate flcHuber(quadratic, wLambda1, connectionsSingle, wLambda2Single, 100, accuracy, 100000, lambda1, lambda2, penEnum::Huber, huberParam);
  flcHuber.runAlgorithm();

  FusedLassoCoordinate flc(quadratic, wLambda1, connectionsSingle, wLambda2Single, maxIterInner, accuracy, maxActivateVars, lambda1, lambda2, penEnum::L1);
  flc.runAlgorithm();
  beta = flc.getBetaOriginal(singleFusions);

  getEqualFusions(newFusions, newFusedGroupSize);

  if(areFusionsEqual(fusions, newFusions)) return(false);
  fusions = newFusions;
  fusedGroupSize = newFusedGroupSize;
  fusedCacheValid_ = false;
  return(true);
}


bool FusedLasso::identifyNewFusions(double lastIterChange, FusionStrategy maxFusion) {
  vector<int> newFusions;
  vector<int> newFusedGroupSize;

  if(maxFusion == FusionStrategy::None) return(false);

  const int maxLevel = static_cast<int>(maxFusion);

  if(lastIterChange > accuracy) fusionLevel = 0;
  else ++fusionLevel;

  if(fusionLevel == 0 && fusionLevel <= maxLevel) {
    getEqualFusions(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) ++fusionLevel;
  }
  if(fusionLevel == 1 && fusionLevel <= maxLevel) {
    getSplitFusionsInactive(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) ++fusionLevel;
  }
  if(fusionLevel == 2 && fusionLevel <= maxLevel) {
    getSplitFusionsActive(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) ++fusionLevel;
  }

  if(fusionLevel > maxLevel || fusionLevel == 3) return(false);
  fusions = newFusions;
  fusedGroupSize = newFusedGroupSize;
  fusedCacheValid_ = false;
  return(true);
}


void FusedLasso::rebuildFusedCache_() {
  vector<vector<int>> fusedGroups(fusedGroupSize.size());
  for(size_t i = 0; i < fusions.size(); ++i) fusedGroups[fusions[i]].push_back(i);

  cachedFusedX_ = X.createFusedX(fusedGroups);

  cachedWLambda1Coordinate_.assign(fusedGroupSize.size(), 0);
  for(size_t i = 0; i < fusedGroupSize.size(); ++i)
    for(size_t j = 0; j < fusedGroups[i].size(); ++j)
      cachedWLambda1Coordinate_[i] += wLambda1[fusedGroups[i][j]];

  pg.getFusedConnectionsWeights(fusions, fusedGroupSize.size(), cachedConnectionsFused_, cachedWLambda2Fused_);
  fusedCacheValid_ = true;
}

bool FusedLasso::runFused(penEnum penType) {
  vector<double> betaFused = translateOriginalToFused();

  if(!fusedCacheValid_) rebuildFusedCache_();

  shared_ptr<QuadraticDerivative> quadDer = make_shared<QuadraticDerivativeDiagonal>(cachedFusedX_, y, wObs, betaFused);
  FusedLassoCoordinate flc(quadDer, cachedWLambda1Coordinate_, cachedConnectionsFused_, cachedWLambda2Fused_, maxIterInner, accuracy/10, maxActivateVars, lambda1, lambda2, penType);

  bool result = flc.runAlgorithm();
  innerIterNum += flc.getIterNum();

  beta = flc.getBetaOriginal(fusions);
  for(size_t i = 0; i < beta.size(); ++i) quadratic->updateBeta(i, beta[i]);
  outerIterNum++;
  return result;
}


bool FusedLasso::runFusedGeneral(penEnum penType) {
  if(outerIterNum > maxIterOuter) return(false);

  vector<double> betaFused = translateOriginalToFused();

  if(!fusedCacheValid_) rebuildFusedCache_();

  vector<double> oldBeta;
  double maxBetaChange = accuracy + 1;
  shared_ptr<QuadraticDerivative> quadDer;
  shared_ptr<FusedLassoCoordinate> flc;
  bool result = true;

  while(maxBetaChange > accuracy && result && outerIterNum <= maxIterOuter) {
    oldBeta = betaFused;
    quadDer = make_shared<QuadraticDerivativeLogistic>(cachedFusedX_, y, wObs, betaFused);
    if(quadDer->isExtreme()) return(false);

    flc = make_shared<FusedLassoCoordinate>(quadDer, cachedWLambda1Coordinate_, cachedConnectionsFused_, cachedWLambda2Fused_, maxIterInner, accuracy/10, maxActivateVars, lambda1, lambda2, penType);
    result = flc->runAlgorithm();
    betaFused = flc->getBeta();
    innerIterNum += flc->getIterNum();

    maxBetaChange = maxDiffDoubleVec(oldBeta, betaFused);
    outerIterNum++;
  }

  beta = flc->getBetaOriginal(fusions);
  quadratic = make_shared<QuadraticDerivativeLogistic>(X, y, wObs, beta);
  if(quadDer->isExtreme()) return(false);

  return(result && outerIterNum < maxIterOuter);
}

bool FusedLasso::runAlgorithmL2() {
  outerIterNum = 0;
  innerIterNum = 0;
  fusedGroupSize.resize(p, 1);
  std::iota(fusions.begin(), fusions.end(), 0);

  bool lastRunOK = (regType != regEnum::GAUSSIAN) ? runFusedGeneral(penEnum::L2) : runFused(penEnum::L2);
  return lastRunOK;
}


bool FusedLasso::runAlgorithmHuber() {
  outerIterNum = 0;
  innerIterNum = 0;
  bool lastRunOK;
  vector<double> oldBeta;

  lastRunOK = (regType != regEnum::GAUSSIAN) ? runFusedGeneral(penEnum::L1) : runFused(penEnum::L1);
  if(quadratic->isExtreme()) return false;

  oldBeta = beta;
  identifyNewFusionsHuber();

  vector<int> newFusions = fusions;
  vector<int> newFusedGroupSize = fusedGroupSize;
  double lastIterChange;

  do {
    oldBeta = beta;
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    lastRunOK = (regType != regEnum::GAUSSIAN) ? runFusedGeneral(penEnum::L1) : runFused(penEnum::L1);
    if(quadratic->isExtreme()) return false;
    getEqualFusions(newFusions, newFusedGroupSize);
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    lastIterChange = maxDiffDoubleVec(oldBeta, beta);
  } while(!areFusionsEqual(fusions, newFusions) && lastIterChange > accuracy);

  return(lastRunOK && outerIterNum < maxIterOuter);
}

bool FusedLasso::runAlgorithm(FusionStrategy maxFusion) {
  outerIterNum = 0;
  innerIterNum = 0;
  bool lastRunOK = true;
  bool newFusionDifferent;
  vector<double> oldBeta;
  fusionLevel = 0;

  while(outerIterNum < maxIterOuter && lastRunOK) {
    oldBeta = beta;
    lastRunOK = (regType != regEnum::GAUSSIAN) ? runFusedGeneral(penEnum::L1) : runFused(penEnum::L1);
    if(quadratic->isExtreme()) { lastRunOK = false; break; }

    double lastIterChange = maxDiffDoubleVec(oldBeta, beta);
    newFusionDifferent = identifyNewFusions(lastIterChange, maxFusion);
    if(!newFusionDifferent) break;
  }

  return(lastRunOK && outerIterNum < maxIterOuter);
}

SparseMatrix FusedLasso::runAlgorithm(
    penEnum penType,
    FusionStrategy maxFusion,
    const vector<double>& lambda1Vec,
    const vector<double>& lambda2Vec,
    const int maxNonZero,
    vector<bool>& success,
    vector<int>& outerIterNumVec,
    vector<int>& innerIterNumVec,
    bool verbose) {

  SparseMatrix betaSols(beta.size());
  success.clear(); success.resize(lambda1Vec.size());
  outerIterNumVec.clear(); outerIterNumVec.resize(lambda1Vec.size());
  innerIterNumVec.clear(); innerIterNumVec.resize(lambda1Vec.size());

  for(size_t i = 0; i < lambda1Vec.size(); ++i) {
    setNewLambdas(lambda1Vec[i], lambda2Vec[i]);
    switch(penType) {
    case penEnum::L1:    success[i] = runAlgorithm(maxFusion); break;
    case penEnum::L2:    success[i] = runAlgorithmL2(); break;
    case penEnum::Huber: success[i] = runAlgorithmHuber(); break;
    }
    outerIterNumVec[i] = getOuterIterNum();
    innerIterNumVec[i] = getInnerIterNum();
    betaSols.addColumn(getBeta());
    if(numNonZero(getBeta()) > maxNonZero) break;
    if(!success[i]) break;
  }
  return betaSols;
}

void FusedLasso::setNewLambdas(double lambda1, double lambda2) {
  this->lambda1 = lambda1;
  this->lambda2 = lambda2;
  fusedCacheValid_ = false;
}

double FusedLasso::findMaxLambda1(const vector<int>& exemptVars) {
  vector<vector<int>> connectionsEmpty(p, vector<int>(0));
  vector<vector<double>> wLambda2Empty(p, vector<double>(0));
  vector<double> wLambda1Extreme(p, 1e6);
  for(int v : exemptVars) wLambda1Extreme[v] = 1e-4;

  FusedLassoCoordinate flc(quadratic, wLambda1Extreme, connectionsEmpty, wLambda2Empty, maxIterInner, accuracy, maxActivateVars, 1, 1, penEnum::L1);
  flc.runAlgorithm();
  this->beta = flc.getBeta();

  if(regType != regEnum::GAUSSIAN) {
    quadratic = make_shared<QuadraticDerivativeLogistic>(this->X, this->y, this->wObs, this->beta);
  }

  double highLambda1 = 0;
  for(int i = 0; i < p; ++i) {
    if(quadratic->getBeta(i) == 0) {
      double normDeriv = fabs(quadratic->getDerivative(i) / wLambda1[i]);
      if(normDeriv > highLambda1) highLambda1 = normDeriv;
    }
  }
  return highLambda1;
}


double FusedLasso::findMaxLambda2(double lambda1_here) {
  if(p > 100000) { REprintf("Too big to estimate lambda2; returning 1"); return n; }

  outerIterNum = 0;
  for(int i = 0; i < p; i++) { beta[i] = 0.1; quadratic->updateBeta(i, 0.1); }

  vector<int> newFusions;
  vector<int> newFusedGroupSize;

  lambda1 = lambda1_here;
  lambda2 = n * 2e6;
  setNewLambdas(lambda1, lambda2);
  getEqualFusions(newFusions, newFusedGroupSize, false, 0.001/accuracy);
  runAlgorithmHuber();

  vector<double> oldBeta(p, 0);
  for(int i = 0; i < p; ++i) oldBeta[i] = beta[i];
  bool hasChanged = false;

  lambda2 /= 2;
  for(int i = 0; i < 30; ++i) {
    setNewLambdas(lambda1, lambda2);
    runAlgorithmHuber();
    for(int j = 0; j < p; ++j) {
      if(fabs(oldBeta[j] - beta[j]) > accuracy * 10) hasChanged = true;
    }
    if(hasChanged) return(lambda2 * 2);
    else lambda2 /= 2;
  }
  return(lambda2);
}


void FusedLasso::getBeta(SparseMatrix& betaMat) { betaMat.addColumn(beta); }


void FusedLasso::makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2) {
  pg.makePullAdjustment(beta, nodePull, lambda2, accuracy);
}

double FusedLasso::calcGroupAverage(int pos, vector<double>& nodePull, vector<double>& beta) {
  double curBeta = beta[pos];
  double average = 0;
  double groupSize = 0;
  while(fabs(beta[pos] - curBeta) < accuracy) {
    groupSize++;
    average += nodePull[pos];
    ++pos;
  }
  return average / groupSize;
}


