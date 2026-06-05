/**********************************************************************
***
*** Implements a graph of all the penalties in the fused lasso problems
*** Groups of fused variables will be treated by taking subgraphs
***
***********************************************************************/

#pragma once

#include "FusedLasso_utils.h"


struct GraphEdge {
  GraphEdge* backwards;
  int to;
  double capacity;
  double flow;
};

typedef vector<GraphEdge*> Node;
typedef vector<Node> Nodes;

class Graph
{
public:
  Nodes nodes;
  vector<bool> subNodesICG; // holds a set of nodes for fast querying for function identifyConnectedGroups
  vector<bool> subNodesCWSM; // holds a set of nodes for fast querying for function connectedWithSameValue
  vector<bool> subNodesCT; // for querying in function connectedTo
  vector<bool> subNodes; // holds the set of active nodes for maxFlow
  vector<int> maxFlowNodeList;

  // objects needed for max-flow calculations
  vector<double> sourceCapacity;
  vector<double> sourceFlow;
  vector<double> sinkCapacity;
  vector<double> sinkFlow;

  vector<double> exFlow; // saves the excess flow
  vector<int> dist; // distance label from the sink
  vector<vector<int>> activeByDist; // active nodes by distance, used as LIFO stacks per level
  int level; // largest distance of an active node

  // ========== FUNCTIONS =======================

  // Helper Functions for maximum flow

  // function performs a breadth first search; return value is the distance of each node
  // if from is true, then it will calculate the distance from the start, otherwise the distance to start
  void distance(const bool from=false);

  // push pushes between regular nodes
  // the other push to the source and sink
  bool push(const int from, GraphEdge* e);
  bool pushToSource(const int from);
  bool pushToSink(const int from);

  // finds the new label for node nodeNum
  int findDist(int nodeNum);
  // preprocess the graph; sets the distances from sink and pushes maximum flow from the source
  void preprocess();
  // push/relabel procedure; return value indicates if node is still active
  bool pushRelabel(const int i);
  // return the next largest active node; will return false if there are no more; node will be saved in nodeNum
  bool getLargestActiveNode(int& nodeNum);
  // insert an active node into the data structure
  void insertActiveNode(const int nodeNum);

  // Helper functions for graph

  // adds an edge to the graph (in both nodes); only intended for use by constructors
  void addEdge(const int from, const int to, const double capacityFromTo, const double capacityToFrom, const double flowFromTo = 0, const double flowToFrom = 0);

  // given a set of nodes, return the set of nodes these nodes are connected to
  // excluding the nodes in the input set itself; these are the nodes a subgraph is connected to
  vector<int> connectedTo(const vector<int>& nodeList);

  void breadthFirstSearch(const int startNode, vector<int>& parentTable, vector<int>& parentEdgeIndex, vector<int>& nodeWithParent, int& endNode, vector<double>& maxFlow);

  double addFlow(const int startNode, const int endNode, const vector<int>& parentTable, const vector<int>& parentEdgeIndex, const vector<double>& maxFlow);

  void constructorWorker(const vector<vector<int>>& conn, const vector<vector<double>>& weights);

  // constructor that uses an R object to build the graph
  Graph(Rcpp::List connList, Rcpp::List connWeights);
  Graph(){};

  // destructor that frees all the edges
  ~Graph();

  // add an unconnected node (used for the intercept)
  void addNode();

  // copy constructor
  void makeCopy(const Graph& pg);
  Graph(const Graph& pg);
  Graph& operator= (const Graph& pg);

  // given a set of nodes, divide them into groups that are connected in the current graph
  vector<vector<int>> identifyConnectedGroups(const vector<int>& x);

  // given a node and a vector with values for every node, find the
  // other nodes connected to it that also have the same value
  vector<int> connectedWithSameValue(int node, const vector<double>& values, double accuracy);

  // for the given group, use the pulls and the penalty graph to divide up into smaller groups
  // items in nodePull will stay the same after return
  vector<vector<int>> splitGroup(const vector<int>& groupNodes, vector<double>& nodePull, double adjustment, double lambda2, bool invert);

  // the nodes are returned in external notation
  vector<int> reachableFromSource();
  // get the complement of the set x w.r.t. y
  vector<int> getComplement(vector<int>& x, vector<int>& y);

  // returns true if the algorithm finished successfully
  bool findMaxFlowPrelabelPush();
  bool findMaxFlowEdmondsKarp();

  // sets all the flows to 0
  void initializeMaxFlow(const vector<int>& groupNodes, const vector<double>& nodePull);
  // count the flow going out of the source
  double getFlowFromSource();
  // count the flow going into the sink
  double getFlowIntoSink();

  // get a set with all the nodes
  vector<int> allNodes();

  // get the connections and weights of the fused graph; parameters get overwritten
  void getFusedConnectionsWeights(const vector<int>& fusions, int numFusions, vector<vector<int>>& connections, vector<vector<double>>& weights);

  void makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2, double accuracy);

  // Helper: Initialize all sub-nodes vectors to false
  void initializeSubNodes();

  // Helper: Clear sub-nodes vectors for a given node set
  void clearSubNodesICG(const vector<int>& nodeList);
  void clearSubNodesCWSM(const vector<int>& nodeList);
  void clearSubNodesCT(const vector<int>& nodeList);

  // prints out the whole graph; used for troubleshooting
  void printGraph(ostream& outStream) const;

  void printIntList(ostream& outStream, const vector<int>& x);
};
