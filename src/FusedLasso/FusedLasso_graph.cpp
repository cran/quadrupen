#include "FusedLasso_graph.h"

using namespace std;

void Graph::distance(const bool fromSource)
{
  // set the distances to maximum
  for(int node : maxFlowNodeList) {
    dist[node] = (int)(maxFlowNodeList.size() + 2);
  }
  queue<int> next;
  int u;
  GraphEdge* e;
  Node::iterator edgeIt;

  for(int node : maxFlowNodeList) {
    if(fromSource) {
      if(sourceCapacity[node] - tolerance > sourceFlow[node]) {
        dist[node] = 1;
        next.push(node);
      }
    }
    else {
      if(sinkCapacity[node] - tolerance > sinkFlow[node]) {
        dist[node] = 1;
        next.push(node);
      }
    }
  }

  while(!next.empty())
  {
    u = next.front();
    next.pop();
    for(edgeIt = nodes[u].begin(); edgeIt != nodes[u].end(); ++edgeIt)
    {
      if(!subNodes[(*edgeIt)->to]) continue;
      if(fromSource) e = *edgeIt;
      else e = (*edgeIt)->backwards;
      if(e->flow < e->capacity - tolerance)
      {
        if(dist[(*edgeIt)->to] > dist[u] + 1)
        {
          dist[(*edgeIt)->to] = dist[u] + 1;
          next.push((*edgeIt)->to);
        }
      }
    }
  }
}

bool Graph::push(const int from, GraphEdge *e)
{
  bool isToActive;
  double addFlowBack = min(exFlow[from], e->backwards->flow);
  double addFlow = min(exFlow[from]-addFlowBack, e->capacity - e->flow);
  e->flow += addFlow;
  e->backwards->flow -= addFlowBack;

  exFlow[from] -= addFlow + addFlowBack;
  isToActive = (exFlow[e->to] > tolerance);
  exFlow[e->to] += addFlow + addFlowBack;
  if(!isToActive) insertActiveNode(e->to);
  return(exFlow[from] > tolerance);
}


bool Graph::pushToSource(const int from) {
  double addFlowBack = min(exFlow[from], sourceFlow[from]);
  exFlow[from] -= addFlowBack;
  sourceFlow[from] -= addFlowBack;
  return(exFlow[from] > tolerance);
}

bool Graph::pushToSink(const int from) {
  double addFlow = min(exFlow[from], sinkCapacity[from] - sinkFlow[from]);
  exFlow[from] -= addFlow;
  sinkFlow[from] += addFlow;
  return(exFlow[from] > tolerance);
}


int Graph::findDist(int nodeNum)
{
  Node::iterator edgeIt;

  int newDist;
  if(sourceFlow[nodeNum] > tolerance) {
    newDist = (int)(maxFlowNodeList.size() + 2);
  }
  else {
    newDist = (int)(2 * maxFlowNodeList.size() + 2);
  }

  for(edgeIt = nodes[nodeNum].begin(); edgeIt != nodes[nodeNum].end(); ++edgeIt)
  {
    if(!subNodes[(*edgeIt)->to]) continue;
    if((*edgeIt)->backwards->flow + (*edgeIt)->capacity - (*edgeIt)->flow > tolerance)
    {
      newDist = min(newDist, dist[(*edgeIt)->to]+1);
    }
  }
  return(newDist);
}


void Graph::preprocess()
{
  activeByDist.assign(2 * maxFlowNodeList.size() + 3, vector<int>());
  level = -1;

  distance();

  for(int node : maxFlowNodeList)
  {
    exFlow[node] = sourceCapacity[node] - sourceFlow[node];
    sourceFlow[node] = sourceCapacity[node];
    if(exFlow[node] > tolerance) insertActiveNode(node);
  }
}

bool Graph::pushRelabel(const int i)
{
  bool didPush = false;
  Node::iterator edgeIt;

  if(sinkCapacity[i] - sinkFlow[i] > tolerance) {
    if(!pushToSink(i)) return(false);
    didPush = true;
  }

  if(dist[i] > (int)(maxFlowNodeList.size() + 1) && sourceFlow[i] > 0) {
    if(!pushToSource(i)) return(false);
    didPush = true;
  }

  for(edgeIt = nodes[i].begin(); edgeIt != nodes[i].end(); ++edgeIt)
  {
    if(!subNodes[(*edgeIt)->to]) continue;
    if((dist[i] > dist[(*edgeIt)->to]) && ((*edgeIt)->capacity > (*edgeIt)->flow + tolerance))
    {
      didPush = true;
      if(!push(i, *edgeIt)) return(false);
    }
  }
  if(!didPush) dist[i] = findDist(i);

  return(true);
}


bool Graph::getLargestActiveNode(int& nodeNum)
{
  if(level < 0) return(false);
  if(activeByDist[level].empty())
  {
    do { --level; }
    while((level >= 0) && (activeByDist[level].empty()));
    if(level < 0) return(false);
  }
  nodeNum = activeByDist[level].back();
  activeByDist[level].pop_back();
  return(true);
}

void Graph::insertActiveNode(const int nodeNum)
{
  if(dist[nodeNum] > level) level = dist[nodeNum];
  activeByDist[dist[nodeNum]].push_back(nodeNum);
}


bool Graph::findMaxFlowPrelabelPush()
{
  int activeNodeNum;
  preprocess();
  while(getLargestActiveNode(activeNodeNum))
  {
    if(pushRelabel(activeNodeNum)) insertActiveNode(activeNodeNum);
  }
  return(true);
}


bool Graph::findMaxFlowEdmondsKarp() {
  vector<int> parentTable(nodes.size(), -1);
  vector<int> parentEdgeIndex(nodes.size(), -1);
  vector<int> nodeWithParent;
  vector<double> maxFlow(nodes.size(), 0);
  int endNode;
  double lastFlow;

  for(int node : maxFlowNodeList) {
    lastFlow = 1;
    while(lastFlow > 0) {
      breadthFirstSearch(node, parentTable, parentEdgeIndex, nodeWithParent, endNode, maxFlow);
      lastFlow = addFlow(node, endNode, parentTable, parentEdgeIndex, maxFlow);
    }
  }
  return(true);
}


void Graph::addEdge(const int from, const int to, const double capacityFromTo, const double capacityToFrom, const double flowFromTo, const double flowToFrom)
{
  GraphEdge* e1 = new(GraphEdge);
  GraphEdge* e2 = new(GraphEdge);

  e1->to = to;   e2->to = from;
  e1->backwards = e2; e2->backwards = e1;
  e1->flow = flowFromTo; e2->flow = flowToFrom;
  e1->capacity = capacityFromTo; e2->capacity = capacityToFrom;

  nodes[from].push_back(e1);
  nodes[to].push_back(e2);
}


vector<int> Graph::connectedTo(const vector<int>& nodeList)
{
  vector<int> conn;
  Node::iterator edgeIt;

  for(int node : nodeList) subNodesCT[node] = true;

  for(int node : nodeList)
  {
    if((size_t)node < nodes.size()) {
      for(edgeIt = nodes[node].begin(); edgeIt != nodes[node].end(); ++edgeIt) {
        if(!subNodesCT[(*edgeIt)->to]) conn.push_back((*edgeIt)->to);
      }
    }
  }

  clearSubNodesCT(nodeList);

  sort(conn.begin(), conn.end());
  conn.erase(unique(conn.begin(), conn.end()), conn.end());
  return(conn);
}


void Graph::breadthFirstSearch(const int startNode, vector<int>& parentTable, vector<int>& parentEdgeIndex, vector<int>& nodeWithParent, int& endNode, vector<double>& maxFlow) {
  vector<int> nodesToSearch;

  for(int n : nodeWithParent) parentTable[n] = -1;
  nodeWithParent.clear();

  nodesToSearch.push_back(startNode);
  maxFlow[startNode] = sourceCapacity[startNode] - sourceFlow[startNode];

  if(maxFlow[startNode] == 0) { endNode = startNode; return; }

  int u;
  GraphEdge* e;
  while(nodesToSearch.size() > 0) {
    u = nodesToSearch.back();
    if(sinkCapacity[u] - sinkFlow[u] > 0) { endNode = u; return; }

    nodesToSearch.pop_back();
    for(size_t edgeIndex = 0; edgeIndex < nodes[u].size(); ++edgeIndex) {
      e = nodes[u][edgeIndex];
      if(!subNodes[e->to]) continue;
      double addFlow = e->capacity - e->flow + e->backwards->flow;
      if(addFlow > 0 && parentTable[e->to] == -1) {
        parentTable[e->to] = u;
        parentEdgeIndex[e->to] = edgeIndex;
        nodeWithParent.push_back(e->to);
        maxFlow[e->to] = min(maxFlow[u], addFlow);
        nodesToSearch.push_back(e->to);
      }
    }
  }

  endNode = startNode;
  maxFlow[startNode] = 0;
}


double Graph::addFlow(const int startNode, const int endNode, const vector<int>& parentTable, const vector<int>& parentEdgeIndex, const vector<double>& maxFlow) {
  double totalFlow = min(maxFlow[endNode], sinkCapacity[endNode] - sinkFlow[endNode]);

  sourceFlow[startNode] += totalFlow;
  sinkFlow[endNode] += totalFlow;
  if(startNode == endNode) return totalFlow;

  int curNode = endNode;
  int parentEdge, parentNode;
  GraphEdge* e;
  while(curNode != startNode) {
    parentNode = parentTable[curNode];
    parentEdge = parentEdgeIndex[curNode];
    e = nodes[parentNode][parentEdge];

    double restFlow = totalFlow;
    if(e->backwards->flow > 0) {
      if(e->backwards->flow < totalFlow) {
        restFlow -= e->backwards->flow;
        e->backwards->flow = 0;
        e->flow += restFlow;
      }
      else { e->backwards->flow -= totalFlow; }
    }
    else { e->flow += totalFlow; }

    curNode = parentNode;
  }

  return totalFlow;
}


void Graph::printGraph(ostream& outStream) const
{
  Node::const_iterator edgeIt;

  for(int node : maxFlowNodeList) {
    if(sourceCapacity[node] > 0)
      outStream << node << ": Source Cap: " << sourceCapacity[node] << " Flow: " << sourceFlow[node] << endl;
    if(sinkCapacity[node] > 0)
      outStream << node << ": Sink Cap: " << sinkCapacity[node] << " Flow: " << sinkFlow[node] << endl;
  }

  for(size_t node = 0; node < nodes.size(); ++node)
  {
    outStream << "Node Number: " << node << endl;
    outStream << "Dist: " << dist[node] << endl;
    outStream << "ExFlow: " << exFlow[node] << endl;
    outStream << "Edges:" << endl;
    for(edgeIt = nodes[node].begin(); edgeIt != nodes[node].end(); ++edgeIt)
      outStream << "To: " << (*edgeIt)->to << " Cap: " << (*edgeIt)->capacity << " Flow: " << (*edgeIt)->flow << endl;
    outStream << endl;
  }
  outStream << endl;
}


void Graph::constructorWorker(const vector<vector<int>>& conn, const vector<vector<double>>& weights) {
  size_t numberOfNodes = conn.size();

  nodes.clear(); nodes.resize(numberOfNodes);
  sourceCapacity.assign(numberOfNodes, 0);
  sourceFlow.assign(numberOfNodes, 0);
  sinkCapacity.assign(numberOfNodes, 0);
  sinkFlow.assign(numberOfNodes, 0);
  subNodesICG.assign(numberOfNodes, false);
  subNodesCWSM.assign(numberOfNodes, false);
  subNodesCT.assign(numberOfNodes, false);
  subNodes.assign(numberOfNodes, false);
  exFlow.assign(numberOfNodes, 0);
  dist.assign(numberOfNodes, 0);
  activeByDist.clear();
  level = 0;
  maxFlowNodeList.clear();

  for(size_t i = 0; i < numberOfNodes; ++i)
  {
    for(size_t j = 0; j < conn[i].size(); ++j)
    {
      int node2 = conn[i][j];
      if(node2 > (int)i) addEdge(i, node2, weights[i][j], weights[i][j]);
    }
  }
}


Graph::Graph(Rcpp::List connList, Rcpp::List connWeights)
{
  int numberOfNodes = connList.length();

  vector<vector<int>> conn(numberOfNodes);
  vector<vector<double>> weights(numberOfNodes);

  for(int i = 0; i < numberOfNodes; ++i)
  {
    Rcpp::IntegerVector connOneNode = connList[i];
    Rcpp::NumericVector weightOneNode = connWeights[i];
    int numOfConn = connOneNode.length();
    conn[i].resize(numOfConn);
    weights[i].resize(numOfConn);
    for(int j = 0; j < numOfConn; ++j) {
      conn[i][j] = connOneNode[j] - 1;
      weights[i][j] = weightOneNode[j];
    }
  }

  constructorWorker(conn, weights);
}


Graph::~Graph()
{
  Nodes::iterator nodeIt;
  Node::iterator edgeIt;
  for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt) {
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) {
      delete(*edgeIt);
      *edgeIt = NULL;
    }
  }
}

void Graph::addNode() {
  nodes.push_back(Node(0));
  subNodesICG.push_back(false);
  subNodesCWSM.push_back(false);
  subNodesCT.push_back(false);
  subNodes.push_back(false);
  sourceCapacity.push_back(0);
  sinkCapacity.push_back(0);
  sourceFlow.push_back(0);
  sinkFlow.push_back(0);
  exFlow.push_back(0);
  dist.push_back(0);
}


void Graph::makeCopy(const Graph& pg) {
  Node::const_iterator edgeIt;
  GraphEdge* edgePtr;
  GraphEdge* edgePtrBack;

  subNodesICG = pg.subNodesICG;
  subNodesCWSM = pg.subNodesCWSM;
  subNodes = pg.subNodes;
  subNodesCT = pg.subNodesCT;
  maxFlowNodeList = pg.maxFlowNodeList;
  sourceCapacity = pg.sourceCapacity;
  sourceFlow = pg.sourceFlow;
  sinkCapacity = pg.sinkCapacity;
  sinkFlow = pg.sinkFlow;
  exFlow = pg.exFlow;
  dist = pg.dist;
  activeByDist = pg.activeByDist;
  level = pg.level;

  nodes.clear(); nodes.resize(pg.nodes.size(), vector<GraphEdge*>());

  for(size_t node = 0; node < pg.nodes.size(); ++node) {
    for(edgeIt = pg.nodes[node].begin(); edgeIt != pg.nodes[node].end(); ++edgeIt) {
      if((size_t)(*edgeIt)->to > node) {
        edgePtr = *edgeIt;
        edgePtrBack = edgePtr->backwards;
        addEdge(node, edgePtr->to, edgePtr->capacity, edgePtrBack->capacity, edgePtr->flow, edgePtrBack->flow);
      }
    }
  }
}


Graph::Graph(const Graph& pg) { makeCopy(pg); }

Graph& Graph::operator= (const Graph& pg) { makeCopy(pg); return(*this); }


vector<vector<int>> Graph::identifyConnectedGroups(const vector<int>& x) {
  vector<int> currentGroup;
  vector<int> stack;
  vector<vector<int>> allGroups;

  for(int node : x) subNodesICG[node] = true;

  for(int startNode : x) {
    if(subNodesICG[startNode]) {
      subNodesICG[startNode] = false;
      currentGroup.clear();
      currentGroup.push_back(startNode);
      stack.push_back(startNode);
      while(!stack.empty()) {
        vector<int> neighbours = connectedTo(stack);
        stack.clear();
        for(int nb : neighbours) {
          if(subNodesICG[nb]) {
            subNodesICG[nb] = false;
            currentGroup.push_back(nb);
            stack.push_back(nb);
          }
        }
      }
      allGroups.push_back(currentGroup);
    }
  }
  clearSubNodesICG(x);
  return(allGroups);
}


vector<int> Graph::connectedWithSameValue(int node, const vector<double>& values, double accuracy) {
  vector<int> stack;
  vector<int> connected;
  double curVal = values[node];

  stack.push_back(node);
  connected.push_back(node);
  subNodesCWSM[node] = true;

  while(!stack.empty()) {
    vector<int> stackConnect = connectedTo(stack);
    stack.clear();
    for(int nb : stackConnect) {
      bool valuesEqual = (curVal == 0 && values[nb] == 0) ||
                         (curVal != 0 && fabs(values[nb] - curVal) < accuracy/10);
      if(valuesEqual && !subNodesCWSM[nb]) {
        connected.push_back(nb);
        subNodesCWSM[nb] = true;
        stack.push_back(nb);
      }
    }
  }

  clearSubNodesCWSM(connected);
  return connected;
}


vector<vector<int>> Graph::splitGroup(const vector<int>& groupNodes, vector<double>& nodePull, double adjustment, double lambda2, bool invert) {
  // Save nodePull values, then transform in-place for max-flow
  vector<double> savedPull(groupNodes.size());
  for(size_t i = 0; i < groupNodes.size(); ++i) {
    savedPull[i] = nodePull[groupNodes[i]];
    nodePull[groupNodes[i]] += adjustment;
    if(invert) nodePull[groupNodes[i]] *= -1;
    nodePull[groupNodes[i]] /= lambda2;
  }

  initializeMaxFlow(groupNodes, nodePull);
  findMaxFlowEdmondsKarp();
  vector<int> mfgNodes = reachableFromSource();
  vector<vector<int>> foo = identifyConnectedGroups(mfgNodes);

  for(size_t i = 0; i < groupNodes.size(); ++i) {
    nodePull[groupNodes[i]] = savedPull[i];
  }
  return foo;
}

vector<int> Graph::reachableFromSource()
{
  vector<int> reachable;
  distance(true);
  for(int node : maxFlowNodeList) {
    if(dist[node] <= (int)maxFlowNodeList.size()) reachable.push_back(node);
  }
  return(reachable);
}


vector<int> Graph::getComplement(vector<int>& x, vector<int>& y)
{
  vector<int> complement;
  sort(x.begin(), x.end());
  sort(y.begin(), y.end());

  auto xIt = x.begin();
  for(int yVal : y) {
    while(xIt != x.end() && *xIt < yVal) ++xIt;
    if(xIt == x.end() || yVal != *xIt) complement.push_back(yVal);
  }
  return(complement);
}


void Graph::initializeMaxFlow(const vector<int>& groupNodes, const vector<double>& nodePull) {
  Node::iterator edgeIt;

  for(int node : maxFlowNodeList) subNodes[node] = false;

  maxFlowNodeList = groupNodes;

  for(int node : groupNodes) {
    for(edgeIt = nodes[node].begin(); edgeIt != nodes[node].end(); ++edgeIt)
      (*edgeIt)->flow = 0;
    sourceFlow[node] = 0;
    sinkFlow[node] = 0;
    exFlow[node] = 0;
    dist[node] = 0;
    subNodes[node] = true;
    if(nodePull[node] > 0) {
      sourceCapacity[node] = nodePull[node];
      sinkCapacity[node] = 0;
    }
    else {
      sourceCapacity[node] = 0;
      sinkCapacity[node] = -nodePull[node];
    }
  }

  activeByDist.clear();
  level = -1;
}


double Graph::getFlowFromSource() {
  double res = 0;
  for(int node : maxFlowNodeList) res += sourceFlow[node];
  return res;
}

double Graph::getFlowIntoSink() {
  double res = 0;
  for(int node : maxFlowNodeList) res += sinkFlow[node];
  return res;
}


vector<int> Graph::allNodes()
{
  vector<int> all;
  all.reserve(nodes.size());
  for(size_t i = 0; i < nodes.size(); ++i) all.push_back((int)i);
  return(all);
}


void Graph::getFusedConnectionsWeights(const vector<int>& fusions, int numFusions, vector<vector<int>>& connections, vector<vector<double>>& weights) {
  Node::iterator nodeIt;

  connections.clear(); connections.resize(numFusions);
  weights.clear(); weights.resize(numFusions);

  vector<vector<int>> mapFusedToOriginal(numFusions);
  for(size_t i = 0; i < fusions.size(); ++i) {
    mapFusedToOriginal[fusions[i]].push_back((int)i);
  }

  vector<int> usedFusedGroups;
  vector<int> fusedGroupPosMap(numFusions, -1);
  int toGroup;

  for(int curGrp = 0; curGrp < numFusions; ++curGrp) {
    for(int origNode : mapFusedToOriginal[curGrp]) {
      for(nodeIt = nodes[origNode].begin(); nodeIt != nodes[origNode].end(); ++nodeIt) {
        toGroup = fusions[(*nodeIt)->to];
        if(toGroup != curGrp) {
          if(fusedGroupPosMap[toGroup] == -1) {
            fusedGroupPosMap[toGroup] = usedFusedGroups.size();
            usedFusedGroups.push_back(toGroup);
            connections[curGrp].push_back(toGroup);
            weights[curGrp].push_back(0);
          }
          weights[curGrp][fusedGroupPosMap[toGroup]] += (*nodeIt)->capacity;
        }
      }
    }
    for(int g : usedFusedGroups) fusedGroupPosMap[g] = -1;
    usedFusedGroups.clear();
  }
}


void Graph::makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2, double accuracy) {
  Node::iterator nodeIt;
  for(int pos = 0; pos < (int)nodePull.size(); ++pos) {
    for(nodeIt = nodes[pos].begin(); nodeIt != nodes[pos].end(); ++nodeIt) {
      if(beta[pos] - beta[(*nodeIt)->to] > accuracy)
        nodePull[pos] += lambda2 * (*nodeIt)->capacity;
      else if(beta[pos] - beta[(*nodeIt)->to] < -accuracy)
        nodePull[pos] -= lambda2 * (*nodeIt)->capacity;
    }
  }
}


void Graph::printIntList(ostream& outStream, const vector<int>& x) {
  for(int val : x) outStream << val << ", ";
  outStream << endl;
}


// ============================================================================
// HELPER FUNCTIONS FOR SUBNODE MANAGEMENT
//

void Graph::initializeSubNodes() {
  int size = nodes.size();
  subNodesICG.assign(size, false);
  subNodesCWSM.assign(size, false);
  subNodesCT.assign(size, false);
  subNodes.assign(size, false);
}

void Graph::clearSubNodesICG(const vector<int>& nodeList) {
  for(int node : nodeList) subNodesICG[node] = false;
}

void Graph::clearSubNodesCWSM(const vector<int>& nodeList) {
  for(int node : nodeList) subNodesCWSM[node] = false;
}

void Graph::clearSubNodesCT(const vector<int>& nodeList) {
  for(int node : nodeList) subNodesCT[node] = false;
}
