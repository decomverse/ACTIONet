#ifndef DAG
#define DAG

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "util.h"

using namespace std;

class DAGNode {
 public:
  inline DAGNode(){};

  inline DAGNode(string node_name, unsigned int node_id,
                 map<string, unsigned>& geneNamesToIDs) {
    name = node_name;
    id = node_id;
    descendents.insert(node_id);
    addDescendentToVec(node_id);

    if (geneNamesToIDs.count(node_name) == 1) {
      genes.insert(geneNamesToIDs[node_name]);
      isTermGene = true;
      weight = 0;
    } else {
      isTermGene = false;
    }

    // if (node_name.substr(0,3) != "S00") {
    /*if ((node_name.substr(0,3) != "S00") && (node_name.substr(0,5) != "GRMZM")
    && (node_name.substr(0,2) != "AC")
        && (node_name.substr(0,2) != "EF") && (node_name.substr(0,2) != "AF") &&
    (node_name.substr(0,4) != "TCGA")) { isTermGene = false; } else { isTermGene
    = true; weight = 0; if (geneNamesToIDs.count(node_name) == 1) {
        genes.insert(geneNamesToIDs[node_name]);
      } else {
        unsigned newGeneID = geneNamesToIDs.size();
        geneNamesToIDs[node_name] = newGeneID;
        genes.insert(newGeneID);
      }
      }*/
  }

  inline DAGNode(string node_name, unsigned int node_id) {
    name = node_name;
    id = node_id;
    descendents.insert(node_id);
    addDescendentToVec(node_id);
    if (node_name.substr(0, 3) != "S00") {
      isTermGene = false;
    } else {
      isTermGene = true;
      weight = 0;
    }
  }

  inline void addParent(unsigned int newParent) { parents.insert(newParent); }

  inline set<unsigned int>::iterator getParentsBegin() {
    return parents.begin();
  }

  inline set<unsigned int>::iterator getParentsEnd() { return parents.end(); }

  inline const set<unsigned int>& getParents() { return parents; }

  inline bool isParent(unsigned possibleParent) {
    return parents.count(possibleParent);
  }

  inline bool isChild(unsigned possibleChild) {
    return children.count(possibleChild);
  }

  inline bool isGeneContained(unsigned possibleGene) {
    return genes.count(possibleGene);
  }

  inline bool isDescendent(unsigned possibleDescendent) {
    // return descendents.count(possibleDescendent);
    if (possibleDescendent < descendentsVec.size()) {
      return descendentsVec[possibleDescendent];
    } else {
      return false;
    }
  }

  inline bool isSibling(unsigned possibleSibling) {
    // return siblings.count(possibleSibling);
    if (possibleSibling < siblingsVec.size()) {
      return siblingsVec[possibleSibling];
    } else {
      return false;
    }
  }

  inline unsigned numParents() { return parents.size(); }

  inline void addChild(unsigned int newChild) { children.insert(newChild); }

  inline set<unsigned int>::iterator getChildrenBegin() {
    return children.begin();
  }

  inline set<unsigned int>::iterator getChildrenEnd() { return children.end(); }

  inline const set<unsigned int>& getChildren() { return children; }

  inline unsigned numChildren() { return children.size(); }

  // This only inserts genes into node's gene list.  To maintain gene lists of
  // all genes in graph, use DAGraph::addGenesToAncestors
  inline void addGenes(set<unsigned>::iterator genesBegin,
                       set<unsigned>::iterator genesEnd) {
    genes.insert(genesBegin, genesEnd);
  }

  inline set<unsigned>::iterator getGenesBegin() { return genes.begin(); }

  inline set<unsigned>::iterator getGenesEnd() { return genes.end(); }

  inline unsigned numGenes() { return genes.size(); }

  inline void addDescendents(set<unsigned>::iterator descendentsBegin,
                             set<unsigned>::iterator descendentsEnd) {
    for (set<unsigned>::iterator it = descendentsBegin; it != descendentsEnd;
         ++it) {
      addDescendentToVec(*it);
    }
    descendents.insert(descendentsBegin, descendentsEnd);
  }

  inline void addDescendentToVec(unsigned newID) {
    if (newID >= descendentsVec.size()) {
      descendentsVec.resize(newID + 1, false);
    }
    descendentsVec[newID] = true;
  }

  inline void addSibling(unsigned newID) {
    if (newID >= siblingsVec.size()) {
      siblingsVec.resize(newID + 1, false);
    }
    siblingsVec[newID] = true;
    siblings.insert(newID);
  }

  inline set<unsigned>::iterator getSiblingsBegin() { return siblings.begin(); }

  inline set<unsigned>::iterator getSiblingsEnd() { return siblings.end(); }

  inline set<unsigned>::iterator getDescendentsBegin() {
    return descendents.begin();
  }

  inline set<unsigned>::iterator getDescendentsEnd() {
    return descendents.end();
  }

  inline set<unsigned>& getDescendents() { return descendents; }

  inline void addAncestor(unsigned newAncestor) {
    ancestors.insert(newAncestor);
  }

  inline set<unsigned>::iterator getAncestorsBegin() {
    return ancestors.begin();
  }

  inline set<unsigned>::iterator getAncestorsEnd() { return ancestors.end(); }

  inline int numDescendents() { return descendents.size(); }

  inline int numAncestors() { return ancestors.size(); }

  inline string getName() { return name; }

  inline unsigned int getID() { return id; }

  inline bool isGene() { return isTermGene; }

  inline void setWeight(double newWeight) { weight = newWeight; }

  inline double getWeight() { return weight; }

 private:
  string name;
  unsigned int id;
  set<unsigned int> parents;
  set<unsigned int> children;
  set<unsigned int> siblings;
  set<unsigned int> genes;
  set<unsigned int> descendents;
  set<unsigned int> ancestors;
  vector<bool> descendentsVec;
  vector<bool> siblingsVec;
  bool isTermGene;
  double weight;
};

class DAGraph {
 public:
  inline unsigned addNode(string nodeName,
                          map<string, unsigned>& geneNamesToIDs) {
    unsigned int nodeID = nodes.size();
    nodes.push_back(DAGNode(nodeName, nodeID, geneNamesToIDs));
    nodeNamesToIDs[nodeName] = nodeID;
    return nodeID;
  }

  inline unsigned addNode() {
    unsigned int nodeID = nodes.size();
    string nodeName = to_string(nodeID);
    nodes.push_back(DAGNode(nodeName, nodeID));
    nodeNamesToIDs[nodeName] = nodeID;
    return nodeID;
  }

  inline void reserveNodes(unsigned numNodes) {
    // cout << "Reserving " << numNodes << endl;
    nodes.reserve(numNodes);
  }

  inline vector<DAGNode>::iterator nodesBegin() { return nodes.begin(); }

  inline vector<DAGNode>::iterator nodesEnd() { return nodes.end(); }

  inline unsigned numNodes() { return nodes.size(); }

  inline void setWeight(double newWeight, unsigned nodeID) {
    nodes[nodeID].setWeight(newWeight);
  }

  inline double getWeight(unsigned nodeID) {
    if (isGene(nodeID)) {
      return geneWeight;
    }
    return nodes[nodeID].getWeight();
  }

  void addEdge(string parentName, string childName,
               map<string, unsigned>& geneNamesToIDs,
               string edgeType = "default") {
    unsigned int parentID;
    unsigned int childID;
    // cout << parentName << "\t" << childName << endl;
    // If parent or child doesn't already exist, add it
    map<string, unsigned int>::iterator parentIt =
        nodeNamesToIDs.find(parentName);
    if (parentIt == nodeNamesToIDs.end()) {
      addNode(parentName, geneNamesToIDs);
      parentID = getID(parentName);
    } else {
      parentID = parentIt->second;
    }
    map<string, unsigned int>::iterator childIt =
        nodeNamesToIDs.find(childName);
    if (childIt == nodeNamesToIDs.end()) {
      if ((edgeType == terminalName) &&
          (geneNamesToIDs.count(childName) == 0)) {
        unsigned newGeneID = geneNamesToIDs.size();
        geneNamesToIDs[childName] = newGeneID;
      }
      addNode(childName, geneNamesToIDs);
      childID = getID(childName);
    } else {
      childID = childIt->second;
    }

    addEdge(parentID, childID, edgeType);
  }

  // This version does not include error checking to ensure both nodes already
  // exist. It should only be used to add edges between existing nodes.
  inline void addEdge(unsigned parentID, unsigned childID,
                      string edgeType = "default") {
    // Update sibling lists
    /*for (set<unsigned>::iterator previousChildrenIt =
      getChildrenBegin(parentID); previousChildrenIt !=
      getChildrenEnd(parentID); ++previousChildrenIt) {
      nodes[*previousChildrenIt].addSibling(childID);
      nodes[childID].addSibling(*previousChildrenIt);
      }*/

    nodes[parentID].addChild(childID);
    nodes[childID].addParent(parentID);
    if (isGene(childID)) {
      edgeType = terminalName;
    }
    edgeTypes[make_pair(parentID, childID)] = edgeType;

    // Add child's gene list to parent's gene list and to all of parent's
    // ancestors' gene lists
    addGenesToAncestors(nodes[parentID], nodes[childID].getGenesBegin(),
                        nodes[childID].getGenesEnd());

    // Add child's descendent list to parent's descendent list and to all of
    // parent's ancestors' descendent list
    addDescendentsToAncestors(nodes[parentID],
                              nodes[childID].getDescendentsBegin(),
                              nodes[childID].getDescendentsEnd());
  }

  // Recursively adds gene list to parent and all of its ansectors
  void addGenesToAncestors(DAGNode& parent, set<unsigned>::iterator genesBegin,
                           set<unsigned>::iterator genesEnd) {
    for (set<unsigned int>::iterator grandparentsIt = parent.getParentsBegin();
         grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addGenesToAncestors(nodes[*grandparentsIt], genesBegin, genesEnd);
    }
    parent.addGenes(genesBegin, genesEnd);
  }

  // Recursively adds descendents list to parent and all of its ansectors
  void addDescendentsToAncestors(DAGNode& parent,
                                 set<unsigned>::iterator descendentsBegin,
                                 set<unsigned>::iterator descendentsEnd) {
    for (set<unsigned int>::iterator grandparentsIt = parent.getParentsBegin();
         grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addDescendentsToAncestors(nodes[*grandparentsIt], descendentsBegin,
                                descendentsEnd);
    }
    parent.addDescendents(descendentsBegin, descendentsEnd);
    for (set<unsigned>::iterator descendentIt = descendentsBegin;
         descendentIt != descendentsEnd; ++descendentIt) {
      nodes[*descendentIt].addAncestor(parent.getID());
    }
  }

  /*
  // Recursively adds gene list to parent and all of its ansectors.  Handles
  cycles in graph void addGenesToAncestors(set<unsigned int> & nodesTraversed,
  DAGNode & parent, set<string>::iterator genesBegin, set<string>::iterator
  genesEnd) {

    // If this node has already been traversed, then return
    if (nodesTraversed.count(parent.getID()) == 1) {
      return;
    }

    // Otherwise, add the id to the traversed list and continue
    nodesTraversed.insert(parent.getID());
    for (set<unsigned int>::iterator grandparentsIt = parent.getParentsBegin();
         grandparentsIt != parent.getParentsEnd(); ++grandparentsIt) {
      addGenesToAncestors(nodesTraversed, nodes[*grandparentsIt], genesBegin,
  genesEnd);
    }
    parent.addGenes(genesBegin, genesEnd);
  }
  */

  inline set<unsigned int>::iterator getParentsBegin(unsigned int id) {
    return nodes[id].getParentsBegin();
  }

  inline set<unsigned int>::iterator getParentsEnd(unsigned int id) {
    return nodes[id].getParentsEnd();
  }

  inline set<unsigned int>::iterator getChildrenBegin(unsigned int id) {
    return nodes[id].getChildrenBegin();
  }

  inline set<unsigned int>::iterator getChildrenEnd(unsigned int id) {
    return nodes[id].getChildrenEnd();
  }

  inline set<unsigned>::iterator getGenesBegin(unsigned int id) {
    return nodes[id].getGenesBegin();
  }

  inline set<unsigned>::iterator getGenesEnd(unsigned int id) {
    return nodes[id].getGenesEnd();
  }

  inline bool isParent(unsigned possibleParent, unsigned possibleChild) {
    return nodes[possibleChild].isParent(possibleParent);
  }

  inline bool isChild(unsigned possibleChild, unsigned possibleParent) {
    return nodes[possibleParent].isChild(possibleChild);
  }

  inline bool isGeneContained(unsigned possibleGene, unsigned possibleParent) {
    return nodes[possibleParent].isGeneContained(possibleGene);
  }

  inline set<unsigned>::iterator getSiblingsBegin(unsigned id) {
    return nodes[id].getSiblingsBegin();
  }

  inline set<unsigned>::iterator getSiblingsEnd(unsigned id) {
    return nodes[id].getSiblingsEnd();
  }

  inline bool isDescendent(unsigned possibleDescendent,
                           unsigned possibleAncestor) {
    return nodes[possibleAncestor].isDescendent(possibleDescendent);
  }

  inline set<unsigned>::iterator getDescendentsBegin(unsigned id) {
    return nodes[id].getDescendentsBegin();
  }

  inline set<unsigned>::iterator getDescendentsEnd(unsigned id) {
    return nodes[id].getDescendentsEnd();
  }

  inline set<unsigned>& getDescendents(unsigned id) {
    return nodes[id].getDescendents();
  }

  inline set<unsigned>::iterator getAncestorsBegin(unsigned id) {
    return nodes[id].getAncestorsBegin();
  }

  inline set<unsigned>::iterator getAncestorsEnd(unsigned id) {
    return nodes[id].getAncestorsEnd();
  }

  inline bool areSiblings(unsigned possibleSibling1,
                          unsigned possibleSibling2) {
    return nodes[possibleSibling1].isSibling(possibleSibling2);
  }

  inline string getName(unsigned int id) { return nodes[id].getName(); }

  inline unsigned int getID(string name) { return nodeNamesToIDs[name]; }

  inline bool isGene(unsigned int id) { return nodes[id].isGene(); }

  inline unsigned numGenesInNode(unsigned int id) {
    return nodes[id].numGenes();
  }

  inline DAGNode getNode(unsigned int id) { return nodes[id]; }

  inline string getEdgeType(unsigned parent, unsigned child) {
    return edgeTypes[make_pair(parent, child)];
  }

  inline map<pair<unsigned, unsigned>, string>::iterator edgesBegin() {
    return edgeTypes.begin();
  }

  inline map<pair<unsigned, unsigned>, string>::iterator edgesEnd() {
    return edgeTypes.end();
  }

  inline void setGeneWeight(double newGeneWeight) {
    geneWeight = newGeneWeight;
    return;
  }

  inline void setTerminalName(string newTerminalName) {
    terminalName = newTerminalName;
    return;
  }

  DAGraph() {
    geneWeight = 1;
    terminalName = "gene";
  };

  DAGraph(string fileName, map<string, unsigned>& geneNamesToIDs,
          string terminalNameForThis = "gene") {
    terminalName = terminalNameForThis;
    string line;
    ifstream file(fileName.c_str());
    if (file.is_open()) {
      while (file.good()) {
        getline(file, line);
        vector<string> tokens;
        Utils::Tokenize(line, tokens, "\t");
        if (tokens.size() == 2) {
          addEdge(tokens[0], tokens[1], geneNamesToIDs);
        } else if (tokens.size() >= 3) {
          addEdge(tokens[0], tokens[1], geneNamesToIDs, tokens[2]);
        }
      }
    }
    geneWeight = 1;
  }

 private:
  vector<DAGNode> nodes;
  map<string, unsigned int> nodeNamesToIDs;
  map<pair<unsigned, unsigned>, string> edgeTypes;
  double geneWeight;
  string terminalName;
};

#endif  // DAG
