#ifndef NODE_DISTANCE_OBJECT
#define NODE_DISTANCE_OBJECT

#include <list>
#include "graph_undirected.h"

typedef list<pair< pair< unsigned, unsigned >, double > > sortedDistanceStruct;

bool compDist(pair< pair< unsigned, unsigned >, double > i, pair< pair< unsigned, unsigned >, double > j) { 
  return (i.second > j.second); 
}

class nodeDistanceObject {
 public:

  nodeDistanceObject() { };

  nodeDistanceObject(unsigned numNodes) {
    nodeDistances.reserve(numNodes);
    numNonZeroStartingHere = vector<unsigned>(numNodes, 0);
    totalNonZero = 0;
    for (unsigned i = 0; i < numNodes; ++i) {
      nodeDistances.push_back(vector<double>());
      nodeDistances[i].reserve(numNodes - i);
      for (unsigned j = (i+1); j < numNodes; ++j) {
	double edgeWeight = 0;
	nodeDistances[i].push_back(edgeWeight);
      }
    }
  }

  nodeDistanceObject(graph_undirected & graph) {
    unsigned numNodes = graph.numNodes();
    nodeDistances.reserve(numNodes);
    numNonZeroStartingHere = vector<unsigned>(numNodes, 0);
    totalNonZero = 0;
    for (unsigned i = 0; i < numNodes; ++i) {
      nodeDistances.push_back(vector<double>());
      nodeDistances[i].reserve(numNodes - i);
      for (unsigned j = (i+1); j < numNodes; ++j) {
	double edgeWeight = graph.getEdgeWeight(i,j);
	nodeDistances[i].push_back(edgeWeight);
	if (edgeWeight != 0) {
	  ++numNonZeroStartingHere[i];
	  ++totalNonZero;
	  sortedDistances.push_back(make_pair(make_pair(i,j),edgeWeight));
	}
      }
    }
    sortedDistances.sort(compDist);
  }

  // Function to return subset of original node distance object
  nodeDistanceObject(nodeDistanceObject & original, vector<unsigned> & nodesInSubset) {
    unsigned numNodes = original.numNodes();
    vector<bool> nodeInSubset(false, numNodes);
    for (vector<unsigned>::iterator it = nodesInSubset.begin(); it != nodesInSubset.end(); ++it) {
      nodeInSubset[*it] = true;
    }

    nodeDistances.reserve(numNodes);
    numNonZeroStartingHere = vector<unsigned>(numNodes, 0);
    totalNonZero = 0;
    for (unsigned i = 0; i < numNodes; ++i) {
      nodeDistances.push_back(vector<double>());
      nodeDistances[i].reserve(numNodes - i);
      for (unsigned j = (i+1); j < numNodes; ++j) {
	double edgeWeight = 0;
	if (nodeInSubset[i] && nodeInSubset[j]) {
	  edgeWeight = original.getDistance(nodesInSubset[i],nodesInSubset[j]);
	}
	nodeDistances[i].push_back(edgeWeight);
	if (edgeWeight != 0) {
	  ++numNonZeroStartingHere[i];
	  ++totalNonZero;
	  sortedDistances.push_back(make_pair(make_pair(i,j),edgeWeight));
	}
      }
    }
    sortedDistances.sort(compDist);
  }

  inline double getDistance(unsigned node1, unsigned node2) {
    if (node2 < node1) {
      unsigned tmp = node1;
      node1 = node2;
      node2 = tmp;
    }
    return (nodeDistances[node1][node2-node1-1]);
  }

  inline void setDistance(unsigned node1, unsigned node2, double distance) {
    if (node2 < node1) {
      unsigned tmp = node1;
      node1 = node2;
      node2 = tmp;
    }
    nodeDistances[node1][node2-node1-1] = distance;
    return;
  }

  inline void sortDistances() {
    sortedDistances.clear();
    for (unsigned i = 0; i < numNodes(); ++i) {
      for (unsigned j = (i+1); j < numNodes(); ++j) {
	double edgeWeight = nodeDistances[i][j-i-1];
	if (edgeWeight != 0) {
	  ++numNonZeroStartingHere[i];
	  ++totalNonZero;
	  sortedDistances.push_back(make_pair(make_pair(i,j),edgeWeight));
	}
      }
    }
    sortedDistances.sort(compDist);
  }

  inline unsigned numNodes() {
    return nodeDistances.size();
  }
  
  inline sortedDistanceStruct::iterator sortedDistancesBegin() {
    return sortedDistances.begin();
  }

  inline sortedDistanceStruct::iterator sortedDistancesEnd() {
    return sortedDistances.end();
  }
  
  inline unsigned getTotalNonZero() {
    return totalNonZero;
  }

  inline unsigned getNumNonZeroStartingHere(unsigned i) {
    return numNonZeroStartingHere[i];
  }

 private:
  vector<vector<double> > nodeDistances;
  vector<unsigned> numNonZeroStartingHere;
  unsigned totalNonZero;
  sortedDistanceStruct sortedDistances;
};

#endif //NODE_DISTANCE_OBJECT
