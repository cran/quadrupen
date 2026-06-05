#include "FusedLasso_optim.h"

using namespace std;

FusedLassoCoordinate::FusedLassoCoordinate(shared_ptr<QuadraticDerivative> quadDer, const vector<double>& wLambda1, const vector<vector<int> > &connections, const vector<vector<double> > &wLambda2, int maxIterInner, double accuracy, int maxActivateVars, double lambda1, double lambda2, penEnum penType, double huberParam) : quadratic(quadDer) {
 
    // Set parameters
    this->p = quadDer->getp();
    this->n = quadDer->getn();
    this->lambda1 = lambda1;
    this->lambda2 = lambda2;
    this->accuracy = accuracy;
    this->huberParam = huberParam;
    this->penType = penType;
    this->maxIterInner = maxIterInner;
    this->maxActivateVars = maxActivateVars;
    this->beta = quadDer->getBetaVec();
    
    // Pre-multiply wLambda1 by lambda1
    this->wLambda1Mult.resize(wLambda1.size());
    for(size_t i = 0; i < wLambda1.size(); ++i) {
        wLambda1Mult[i] = wLambda1[i] * lambda1;
    }
    
    // Initialize active variables (all non-zero beta)
    active.clear();
    active.reserve(min(n, p));
    for(int i = 0; i < p; ++i) {
        if(beta[i] != 0) {
            active.push_back(i);
            quadDer->activate(i);
        }
    }
    
    // Store connections and pre-multiply wLambda2 by lambda2
    this->connections = connections;
    this->wLambda2Mult = wLambda2;
    for(size_t i = 0; i < wLambda2Mult.size(); ++i) {
        for(size_t j = 0; j < wLambda2Mult[i].size(); ++j) {
            wLambda2Mult[i][j] *= lambda2;
        }
    }
    
    quadratic = quadDer;
}



bool FusedLassoCoordinate::activateBeta(int pos) {
    // first check if the element is already included
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt == active.end()) {// could not find the element
        active.push_back(pos);
        quadratic->activate(pos);
        return true;
    }
    return false;
}


bool FusedLassoCoordinate::deactivateBeta(int pos) {
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt != active.end()) {
        active.erase(vecIt); // delete the element
        return true;
    }
    return false;
}

void FusedLassoCoordinate::singleStep(int pos) {
    
    // here, view the problem as only having variable pos
    // find the optimal position and move beta there
    double deriv = quadratic->getDerivative(pos);
    double slope = quadratic->getHessian(pos); 

    switch(penType) {
        case penEnum::L1: Steps::findRootL1(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos], rootScratch_); break;
        case penEnum::Huber: Steps::findRootHuber(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos], huberParam, rootScratch_); break;
        case penEnum::L2: Steps::findRootL2(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos]); break;
    }
    quadratic->updateBeta(pos, beta[pos]);
}

double FusedLassoCoordinate::singleIteration(const int iterNum) {
    // Return early if no active variables
    if(active.size() == 0) {
        return 0;
    }

    // Perform single steps for all active variables
    // and track the maximum change
    vector<int>::iterator vecIt;
    double maxNorm = 0;
    
    for(vecIt = active.begin(); vecIt != active.end(); ++vecIt) {
        double oldBeta = beta[*vecIt];
        singleStep(*vecIt);
        
        double change = fabs(beta[*vecIt] - oldBeta);
        if(change > maxNorm) { 
            maxNorm = change; 
        }
    }

    return maxNorm;
}


int FusedLassoCoordinate::activateVariables() {
    vector<double> derivVec = quadratic->getDerivativeVec();
    vector<bool> isActive(beta.size(), false);
    for(int i = 0; i < (int)active.size(); ++i) {
        isActive[active[i]] = true;
    }

    // Collect candidates: (|adjDeriv|, pos) for inactive variables that violate KKT
    vector<pair<double, int>> candidates;
    for(int pos = 0; pos < (int)beta.size(); ++pos) {
        if(!isActive[pos]) {
            double zeroPenalty = wLambda1Mult[pos];
            double adjDeriv = derivVec[pos];
            derivAdjustment(adjDeriv, zeroPenalty, pos);
            if(fabs(adjDeriv) > zeroPenalty) {
                candidates.emplace_back(fabs(adjDeriv), pos);
            }
        }
    }

    // Partial sort: bring the top-maxActivateVars to the front in O(n) average
    int activateCount = (int)std::min((int)candidates.size(), maxActivateVars);
    if(activateCount > 0 && activateCount < (int)candidates.size()) {
        std::nth_element(candidates.begin(), candidates.begin() + activateCount, candidates.end(),
            [](const pair<double,int>& a, const pair<double,int>& b){ return a.first > b.first; });
    }

    for(int i = 0; i < activateCount; ++i) {
        active.push_back(candidates[i].second);
        quadratic->activate(candidates[i].second);
    }

    sort(active.begin(), active.end());
    return activateCount;
}

void FusedLassoCoordinate::derivAdjustment(double& adjDeriv, double& zeroPenalty, int pos) {
    // Adjust derivative based on penalty type and connected variables
    
    switch(penType) {
        case penEnum::L1:
            // L1 penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                double connBeta = beta[connections[pos][i]];
                if(connBeta > 0) {
                    adjDeriv -= wLambda2Mult[pos][i];
                }
                else if(connBeta < 0) {
                    adjDeriv += wLambda2Mult[pos][i];
                }
                else {
                    zeroPenalty += wLambda2Mult[pos][i];
                }
            }
            break;
            
        case penEnum::Huber:
            // Huber penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                double connBeta = beta[connections[pos][i]];
                double threshold = 1.0 / huberParam;
                if(connBeta > threshold) {
                    adjDeriv -= wLambda2Mult[pos][i];
                }
                else if(connBeta < -threshold) {
                    adjDeriv += wLambda2Mult[pos][i];
                }
                else {
                    // Smooth region: quadratic penalty
                    adjDeriv -= huberParam * connBeta * wLambda2Mult[pos][i];
                }
            }
            break;
            
        case penEnum::L2:
            // L2 penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                adjDeriv -= beta[connections[pos][i]] * 2 * wLambda2Mult[pos][i];
            }
            break;
    }
}

int FusedLassoCoordinate::deactivateVariables() {
    // go through all active variables and remove the
    // ones from being active that are equal to 0

    int deleteCount = 0;
    vector<int>::iterator vecIt;
    for(vecIt = active.begin(); vecIt != active.end(); ) {
       if(beta[*vecIt] == 0) { // yes, deactivate
          vecIt = active.erase(vecIt); // due to deletion don't need to move one element ahead
          deleteCount++;
       }
       else {
           vecIt++; 
       }
    }
    return deleteCount;
}

bool FusedLassoCoordinate::runAlgorithm() {
    // Main optimization loop with variable activation
    bool finished = false;
    int newActivations;
    iterNum = 0;
    double curError;

    while(!finished && iterNum < maxIterInner) {
        curError = 1;
        
        // Inner loop: optimize over active variables until convergence
        while(curError > accuracy && iterNum < maxIterInner) {
            curError = singleIteration(iterNum);
            R_CheckUserInterrupt(); // Allow user to cancel computation

            iterNum++;
            if(curError > 1e5) {
                return(false); // Divergence detected
            }
        }

        // Try to activate new variables
        newActivations = activateVariables();

        if(newActivations == 0) {
            // No new activations: algorithm has converged
            finished = true;
        }
        if(iterNum >= maxIterInner) {
            // Maximum iterations reached
            finished = false;
        }
    }

    return(finished);
}

int FusedLassoCoordinate::getIterNum() const {
    return iterNum;
}

vector<double> FusedLassoCoordinate::getBetaOriginal(vector<int>& fusions) {
    vector<double> betaOrig(fusions.size());

    for(unsigned int i = 0; i < fusions.size(); ++i) {
        if((size_t)fusions[i] >= beta.size()) {
            Rcpp::stop("The node has a group that is too large. Node: %d\n", i);
        }
        betaOrig[i] = beta[fusions[i]];
    }
    return betaOrig;
}


