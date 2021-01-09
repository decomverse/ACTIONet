#ifndef GRAPH_UNDIRECTED
#define GRAPH_UNDIRECTED

#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "util.h"

#include <ACTIONet.h>

using namespace std;

inline bool compEdges(pair<pair<unsigned, unsigned>, double> i,
                      pair<pair<unsigned, unsigned>, double> j) {
  return i.second > j.second;
}

inline double logistic_func(double x) { return (1 / (1 + exp(-1 * x))); }

class Node {
 public:
  inline Node(){};

  inline Node(string node_name, unsigned int node_id) {
    name = node_name;
    id = node_id;
  }

  inline void insertInOrder(vector<unsigned>& sortedVector,
                            unsigned newElement) {
    unsigned minPos = 0;
    unsigned maxPos = sortedVector.size() - 1;
    if ((sortedVector.size() == 0) || (newElement > sortedVector[maxPos])) {
      sortedVector.push_back(newElement);
      return;
    }
    if (newElement < sortedVector[0]) {
      sortedVector.insert(sortedVector.begin(), newElement);
      return;
    }
    while ((maxPos - minPos) > 1) {
      unsigned testPos = (minPos + maxPos) / 2;
      if (newElement > sortedVector[testPos]) {
        minPos = testPos;
      } else if (newElement < sortedVector[testPos]) {
        maxPos = testPos;
      }
    }

    vector<unsigned>::iterator it = sortedVector.begin();
    it += maxPos;
    sortedVector.insert(it, newElement);
  }

  inline void addInteractor(unsigned int newInteractor) {
    if (newInteractor >= interactors_bool.size()) {
      interactors_bool.resize(newInteractor + 1, false);
    }
    if (!interactors_bool[newInteractor]) {
      insertInOrder(interactors, newInteractor);
      // interactors.push_back(newInteractor);
      // sort(interactors.begin(), interactors.end());
    }
    interactors_bool[newInteractor] = true;
    // ints_plus_ints_of_ints.insert(newInteractor);
  }

  inline vector<unsigned int>::iterator getInteractorsBegin() {
    return interactors.begin();
  }

  inline vector<unsigned int>::iterator getInteractorsEnd() {
    return interactors.end();
  }

  inline const vector<unsigned int>& getInteractors() { return interactors; }

  inline bool isInteractor(unsigned possibleInteractor) {
    if (possibleInteractor > interactors_bool.size()) {
      return false;
    }
    return interactors_bool[possibleInteractor];
  }

  inline unsigned numInteractors() { return interactors.size(); }

  /*inline set<unsigned int>::iterator getIntsPlusIntsOfIntsBegin() {
    return ints_plus_ints_of_ints.begin();
  }

  inline set<unsigned int>::iterator getIntsPlusIntsOfIntsEnd() {
    return ints_plus_ints_of_ints.end();
    }*/

  inline string getName() { return name; }

  inline unsigned int getID() { return id; }

 private:
  string name;
  unsigned int id;
  vector<unsigned int> interactors;
  vector<bool> interactors_bool;
  // set<unsigned int> ints_plus_ints_of_ints;
};

class graph_undirected {
 public:
  inline int addNode(string nodeName, map<string, unsigned>& nodeNamesToIDs) {
    unsigned int nodeID = nodes.size();
    nodes.push_back(Node(nodeName, nodeID));
    nodeNamesToIDs[nodeName] = nodeID;
    return nodeID;
  }

  inline int addNode(string nodeName, unsigned nodeID) {
    if (nodes.size() <= nodeID) {
      nodes.resize(nodeID + 1);
    }
    nodes[nodeID] = Node(nodeName, nodeID);
    return nodeID;
  }

  inline graph_undirected(unsigned numNodes) {
    keepIntsOfInts = false;
    background = 0;
    nodes.reserve(numNodes);
    for (unsigned i = 0; i < numNodes; ++i) {
      nodes.push_back(Node(to_string(i), i));
    }
  }

  inline vector<Node>::iterator nodesBegin() { return nodes.begin(); }

  inline vector<Node>::iterator nodesEnd() { return nodes.end(); }

  inline unsigned numNodes() { return nodes.size(); }

  void addEdge(string node1Name, string node2Name,
               map<string, unsigned>& nodeNamesToIDs, double edgeWeight = 1) {
    unsigned node1ID;
    unsigned node2ID;

    // If parent or child doesn't already exist, add it
    map<string, unsigned int>::iterator node1It =
        nodeNamesToIDs.find(node1Name);
    if (node1It == nodeNamesToIDs.end()) {
      node1ID = addNode(node1Name, nodeNamesToIDs);
    } else {
      node1ID = node1It->second;
    }
    map<string, unsigned int>::iterator node2It =
        nodeNamesToIDs.find(node2Name);
    if (node2It == nodeNamesToIDs.end()) {
      node2ID = addNode(node2Name, nodeNamesToIDs);
    } else {
      node2ID = node2It->second;
    }

    addEdge(node1ID, node2ID, edgeWeight);
    return;
  }

  // Add edge using node IDs.  Nodes must already exist for this to work
  void addEdge(unsigned node1ID, unsigned node2ID, double edgeWeight = 1) {
    nodes[node1ID].addInteractor(node2ID);
    /*if (keepIntsOfInts) {
      for (set<unsigned>::iterator node2IntsIt = getInteractorsBegin(node2ID);
      node2IntsIt != getInteractorsEnd(node2ID); ++node2IntsIt) {
        nodes[node1ID].addIntOfInt(*node2IntsIt);
        nodes[*node2IntsIt].addIntOfInt(node1ID);
      }
      }*/

    nodes[node2ID].addInteractor(node1ID);
    /*if (keepIntsOfInts) {
      for (set<unsigned>::iterator node1IntsIt = getInteractorsBegin(node1ID);
      node1IntsIt != getInteractorsEnd(node1ID); ++node1IntsIt) {
        nodes[node2ID].addIntOfInt(*node1IntsIt);
        nodes[*node1IntsIt].addIntOfInt(node2ID);
      }
      }*/

    if (node1ID < node2ID) {
      edges[make_pair(node1ID, node2ID)] = edgeWeight;
    } else {
      edges[make_pair(node2ID, node1ID)] = edgeWeight;
    }
    return;
  }

  inline vector<unsigned int>::iterator getInteractorsBegin(unsigned int id) {
    return nodes[id].getInteractorsBegin();
  }

  inline vector<unsigned int>::iterator getInteractorsEnd(unsigned int id) {
    return nodes[id].getInteractorsEnd();
  }

  inline const vector<unsigned int>& getInteractors(unsigned int id) {
    return nodes[id].getInteractors();
  }

  /*inline set<unsigned int>::iterator getIntsPlusIntsOfIntsBegin(unsigned int
  id) { return nodes[id].getIntsPlusIntsOfIntsBegin();
  }

  inline set<unsigned int>::iterator getIntsPlusIntsOfIntsEnd(unsigned int id) {
    return nodes[id].getIntsPlusIntsOfIntsEnd();
    }*/

  inline bool isInteraction(unsigned node1ID, unsigned node2ID) {
    return nodes[node1ID].isInteractor(node2ID);
  }

  inline string getName(unsigned int id) { return nodes[id].getName(); }

  inline unsigned int getID(string name) { return nodeNamesToIDs[name]; }

  inline Node getNode(unsigned int id) { return nodes[id]; }

  inline double getEdgeWeight(unsigned node1, unsigned node2) {
    if (node2 < node1) {
      unsigned tmp = node1;
      node1 = node2;
      node2 = tmp;
    }
    map<pair<unsigned, unsigned>, double>::iterator edgeFinder =
        edges.find(make_pair(node1, node2));
    if (edgeFinder == edges.end()) {
      return background;
    }
    return edgeFinder->second;
  }

  inline double getEdgeWeightLogistic(unsigned node1, unsigned node2) {
    return logistic_func(this->getEdgeWeight(node1, node2));
  }

  inline map<pair<unsigned, unsigned>, double>::iterator edgesBegin() {
    return edges.begin();
  }

  inline map<pair<unsigned, unsigned>, double>::iterator edgesEnd() {
    return edges.end();
  }

  inline unsigned numEdges() { return edges.size(); }

  inline double getBackground() { return background; }

  inline double getBackgroundLogistic() { return logistic_func(background); }

  inline void setBackground(double newBack) { background = newBack; }

  graph_undirected() {
    keepIntsOfInts = false;
    background = 0;
  }

  graph_undirected(map<string, unsigned>& nodeNamesToIDs,
                   bool keepTheseIntsOfInts = false) {
    keepIntsOfInts = keepTheseIntsOfInts;
    nodes.reserve(nodeNamesToIDs.size());
    for (map<string, unsigned>::iterator nodeIt = nodeNamesToIDs.begin();
         nodeIt != nodeNamesToIDs.end(); ++nodeIt) {
      addNode(nodeIt->first, nodeIt->second);
    }
    background = 0;
  }

  graph_undirected(string fileName, map<string, unsigned>& nodeNamesToIDs,
                   bool keepTheseIntsOfInts = false) {
    keepIntsOfInts = keepTheseIntsOfInts;
    nodes.reserve(nodeNamesToIDs.size());
    for (map<string, unsigned>::iterator nodeIt = nodeNamesToIDs.begin();
         nodeIt != nodeNamesToIDs.end(); ++nodeIt) {
      addNode(nodeIt->first, nodeIt->second);
    }

    string line;
    ifstream file(fileName.c_str());
    if (file.is_open()) {
      while (file.good()) {
        getline(file, line);
        vector<string> tokens;
        Utils::Tokenize(line, tokens, "\t");
        if (tokens.size() == 2) {
          addEdge(tokens[0], tokens[1], nodeNamesToIDs);
        } else if (tokens.size() == 3) {
          addEdge(tokens[0], tokens[1], nodeNamesToIDs, stod(tokens[2]));
        }
      }
    }
    background = 0;
  }

  graph_undirected(mat& A) {
    unsigned numNodes = A.n_rows;

    keepIntsOfInts = false;
    background = 0;
    nodes.reserve(numNodes);
    for (unsigned i = 0; i < numNodes; ++i) {
      nodes.push_back(Node(to_string(i), i));
    }

    for (int i = 0; i < numNodes; i++) {
      for (int j = i + 1; j < numNodes; j++) {
        addEdge(i, j, A(i, j));
      }
    }
  }

  inline double maxEdgeWeight() {
    double max = edges.begin()->second;
    for (map<pair<unsigned, unsigned>, double>::iterator edgeIt =
             ++edges.begin();
         edgeIt != edges.end(); ++edgeIt) {
      if (edgeIt->second > max) {
        max = edgeIt->second;
      }
    }
    return max;
  }

  inline double minEdgeWeight() {
    double min = edges.begin()->second;
    for (map<pair<unsigned, unsigned>, double>::iterator edgeIt =
             ++edges.begin();
         edgeIt != edges.end(); ++edgeIt) {
      if (edgeIt->second < min) {
        min = edgeIt->second;
      }
    }
    return min;
  }

  inline double calculateExpectedNumInteractors() {
    double expectedNumInteractors = 0;
    for (vector<Node>::iterator nodeIt = nodes.begin(); nodeIt != nodes.end();
         ++nodeIt) {
      expectedNumInteractors +=
          nodeIt->numInteractors() * nodeIt->numInteractors();
    }
    return expectedNumInteractors / edges.size();
  }

  inline unsigned getNumInteractors(unsigned node) {
    return nodes[node].numInteractors();
  }

  inline double nodeOutDegree(unsigned node) {
    if (node >= nodes.size()) {
      return (nodes.size() * this->getBackground());
    }
    double outDegree = 0;
    for (vector<unsigned>::iterator intsIt = this->getInteractorsBegin(node);
         intsIt != this->getInteractorsEnd(node); ++intsIt) {
      outDegree += this->getEdgeWeight(node, *intsIt);
    }
    unsigned numBackgroundInteractors =
        nodes.size() - getNumInteractors(node) - 1;
    outDegree += numBackgroundInteractors * this->getBackground();
    // cout << nodes[node].getName() << "\t" << outDegree << endl;
    return outDegree;
  }

  inline double nodeOutDegreeLogistic(unsigned node) {
    double outDegree = 0;
    for (vector<unsigned>::iterator intsIt = this->getInteractorsBegin(node);
         intsIt != this->getInteractorsEnd(node); ++intsIt) {
      outDegree += this->getEdgeWeightLogistic(node, *intsIt);
    }
    unsigned numBackgroundInteractors =
        nodes.size() - getNumInteractors(node) - 1;
    outDegree += numBackgroundInteractors * this->getBackgroundLogistic();
    return outDegree;
  }

  inline void netToPowerLaw(double powerLawExponent = 2.0) {
    // Iterate through all edges in the network
    vector<pair<pair<unsigned, unsigned>, double> > edgesVec(edges.begin(),
                                                             edges.end());
    std::sort(edgesVec.begin(), edgesVec.end(), compEdges);
    // unsigned numPossibleEdges = ((nodes.size() * (nodes.size() - 1))/2);
    unsigned numEdges = edgesVec.size();
    double lastVal = edgesVec.begin()->second;
    unsigned numWithThisRank = 0;
    unsigned startOfThisRank = 0;

    for (vector<pair<pair<unsigned, unsigned>, double> >::iterator it =
             edgesVec.begin();
         it != edgesVec.end(); ++it) {
      // cout << getName(it->first.first) << "\t" << getName(it->first.second)
      // << "\t" << it->second << "\t" << lastVal << endl;
      if (it->second == lastVal) {
        ++numWithThisRank;
      } else {
        // cout << "start: " << startOfThisRank << endl;
        // cout << "numWithThisRank: " << numWithThisRank << endl;
        double medianRank = startOfThisRank + (numWithThisRank / 2.0);
        // double transformedScore = pow(1 - (medianRank/numPossibleEdges),
        // powerLawExponent);
        double transformedScore =
            pow(1 - (medianRank / numEdges), powerLawExponent);
        for (unsigned i = 0; i < numWithThisRank; ++i) {
          cout << getName(edgesVec[startOfThisRank + i].first.first) << "\t"
               << getName(edgesVec[startOfThisRank + i].first.second) << "\t"
               << transformedScore << endl;
        }

        // Start counting next rank
        startOfThisRank += numWithThisRank;
        numWithThisRank = 1;
        lastVal = it->second;
      }
    }
    // Last val needs to be considered
    double medianRank = startOfThisRank + (numWithThisRank / 2.0);
    // double transformedScore = pow(1 - (medianRank/numPossibleEdges),
    // powerLawExponent);
    double transformedScore =
        pow(1 - (medianRank / numEdges), powerLawExponent);
    for (unsigned i = 0; i < numWithThisRank; ++i) {
      cout << getName(edgesVec[startOfThisRank + i].first.first) << "\t"
           << getName(edgesVec[startOfThisRank + i].first.second) << "\t"
           << transformedScore << endl;
    }
    return;
  }

 private:
  vector<Node> nodes;
  map<string, unsigned int> nodeNamesToIDs;
  map<pair<unsigned, unsigned>, double> edges;
  bool keepIntsOfInts;
  double background;
};

#endif  // GRAPH_UNDIRECTED
