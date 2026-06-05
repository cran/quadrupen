#include "FusedLasso_debug.h"

using namespace std;

// ── FusedLasso debug methods ─────────────────────────────────────────────────

void FusedLasso::printBetaAndGroup(ostream& out) {
  out << "========== Fusions =================" << endl;
  int curGroup = -1;
  for(int i = 0; i < p; ++i) {
    if(fusions[i] > curGroup) {
      curGroup++;
      out << " Group: " << i << " Beta: " << beta[i];
    }
  }
  out << endl << "========== End Fusions =============" << endl;
}

void FusedLasso::printGroups(vector<vector<int>> x, ostream& out) {
  for(size_t i = 0; i < x.size(); ++i)
    for(int node : x[i]) out << "i: " << i << " : " << node << endl;
}

void FusedLasso::checkSolution(ostream& out) {
  vector<double> nodePull = getPulls();
  vector<double> nodePullOriginal = nodePull;
  vector<double> nodePullAdj = nodePull;
  makePullAdjustment(beta, nodePullAdj, lambda2);
  vector<int> newFusions;
  vector<int> newFusedGroupSize;
  getEqualFusions(newFusions, newFusedGroupSize);

  size_t curPos = 0;
  int curGroup = newFusions[curPos];
  double adjustment;

  out << "Lambda1: " << lambda1 << endl;
  out << "Lambda2: " << lambda2 << endl;

  while(curPos < nodePull.size()) {
    if(beta[curPos] > 0) {
      nodePull[curPos] += lambda1;
      curPos++;
    } else if(beta[curPos] < 0) {
      nodePull[curPos] -= lambda1;
      curPos++;
    } else {
      curGroup = newFusions[curPos];
      adjustment = calcGroupAverage(curPos, nodePullAdj, beta);
      if(fabs(adjustment) > lambda1)
        out << "=============== Problem =================" << endl
            << "Group " << curGroup << " at position " << curPos
            << " has adjustment " << adjustment << endl;
      while(curPos < nodePull.size() && newFusions[curPos] == curGroup) {
        nodePull[curPos] -= adjustment;
        curPos++;
      }
    }
  }

  vector<int> allNodes = pg.allNodes();
  for(size_t i = 0; i < nodePull.size(); ++i) {
    nodePull[i] /= lambda2;
    out << i << ": " << nodePull[i] << endl;
  }

  pg.initializeMaxFlow(allNodes, nodePull);
  pg.findMaxFlowEdmondsKarp();
  GraphEdge* e;
  GraphEdge* eBack;

  for(size_t i = 0; i < beta.size(); ++i) {
    e = pg.nodes[i].back();
    eBack = e->backwards;
    out << "Node " << i << " Beta: " << beta[i] << " Deriv: " << nodePullOriginal[i]
        << " DerivAdj: " << nodePull[i];
    if(i < beta.size()-1) out << "t[i,i+1]: " << e->flow << "Back: " << eBack->flow;
    out << endl;
  }

  out << "================= Whole MaxflowGraph ==================" << endl;
  out << "==================== END WHOLE =======================" << endl;
}

// ── FusedLassoCoordinate debug methods ───────────────────────────────────────

void FusedLassoCoordinate::printBetaActive(ostream& out) {
  out << "===================== Beta Active =====================" << endl;
  for(size_t i = 0; i < beta.size(); ++i)
    if(beta[i] != 0) out << i << ":" << beta[i] << " | ";
  out << endl << "======================== End Beta Active ================" << endl;
}

void FusedLassoCoordinate::printDerivActive(ostream& out) {
  out << "===================== Deriv Active =====================" << endl;
  for(size_t i = 0; i < beta.size(); ++i)
    if(beta[i] != 0) out << i << ":" << quadratic->getDerivative(i) << " | ";
  out << endl << "======================== End Deriv Active ================" << endl;
}

void FusedLassoCoordinate::printBeta(ostream& out) {
  out << "======================= Beta ===============" << endl;
  for(size_t i = 0; i < beta.size(); ++i) out << beta[i] << ", ";
  out << endl << "================ End Beta ================" << endl;
}

void FusedLassoCoordinate::printDerivs(ostream& out) {
  out << "============= Derivs =============" << endl;
  for(size_t i = 0; i < beta.size(); ++i) out << quadratic->getDerivative(i) << " ";
  out << endl << "============= End Derivs =========" << endl;
}

void FusedLassoCoordinate::printConnectionsWeights(ostream& out) {
  out << "============ Conn Weights ============" << endl;
  for(size_t i = 0; i < connections.size(); ++i) {
    out << "Node:" << i << endl;
    for(size_t j = 0; j < connections[i].size(); ++j)
      out << connections[i][j] << " " << wLambda2Mult[i][j] << ";";
    out << endl;
  }
}

void FusedLassoCoordinate::printPosSingleStepInfo(int pos, ostream& out) {
  out << "Info at Step: " << pos << endl;
  out << "Deriv: "    << quadratic->getDerivative(pos) << endl;
  out << "Hessian: "  << quadratic->getHessian(pos) << endl;
  out << "OldBeta: "  << beta[pos] << endl;
  out << "wLambda1: " << wLambda1Mult[pos] << endl;
  out << "Connections: " << endl;
  printVector(connections[pos], out);
  out << "wLambda2Mult: " << endl;
  printVector(wLambda2Mult[pos], out);
}
