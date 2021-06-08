#ifndef DAG_CONSTRUCT
#define DAG_CONSTRUCT

#include <math.h>
#include <time.h>
#include <list>
#include <algorithm>
#include <thread>
#include <mutex>
#include "dag.h"
#include "graph_undirected.h"
#include "graph_undirected_bitset.h"
#include "util.h"
#include "nodeDistanceObject.h"
#include "corrector.h"

void printCluster(const boost::dynamic_bitset<unsigned long> & cluster, vector<string> & nodeIDsToNames) {
	/*
  for (unsigned i = 0; i < cluster.size(); ++i) {
    if (cluster[i]) {
      cout << nodeIDsToNames[i] << ",";
    }
  }
  */
}

bool compPairSecondAscending(const pair< unsigned, unsigned > & i,const pair< unsigned, unsigned > & j) { 
  return (i.second < j.second); 
}


/*bool compNewClustersAndCounts(const pair< unsigned, pair<unsigned, unsigned > > & i,const pair< unsigned, pair<unsigned, unsigned > > & j) {
  if (i.second.first == j.second.first) {
    return (i.second.second < j.second.second);
  }
  return (i.second.first < j.second.first); 
  }*/

bool compNewClustersAndCounts(const pair< unsigned, pair<unsigned, unsigned > > & i,const pair< unsigned, pair<unsigned, unsigned > > & j) {
  if (i.second.second == j.second.second) {
    return (i.second.first > j.second.first);
  }
  return (i.second.second < j.second.second); 
}


bool compPairSecondDescending(const pair< unsigned, unsigned > & i,const pair< unsigned, unsigned > & j) { 
  return (i.second > j.second); 
}

bool compClustersToCombine(const pair<pair<unsigned long, unsigned long>, double> & i, const pair<pair<unsigned long, unsigned long>, double> & j) {
  if (i.second == j.second) {
    if (i.first.first == j.first.first) {
      return i.first.second < j.first.second;
    }
    return i.first.first < j.first.first;
  }
  return i.second > j.second;
}

bool compEdgesToAdd(const pair<pair<unsigned, unsigned>, double> & i, const pair<pair<unsigned, unsigned>, double> & j) {
  if (i.first.first == j.first.first) {
    return i.first.second < j.first.second;
  }
  return i.first.first < j.first.first;
}

bool newClustersComp(const pair<boost::dynamic_bitset<unsigned long>, unsigned> & i, const pair<boost::dynamic_bitset<unsigned long>, unsigned> & j) {
  if (i.second == j.second) {
    return i.first < j.first;
  }
  return i.second < j.second;
}

struct compClusters {
  bool operator() (const vector<unsigned>& i, const vector <unsigned>& j) const {
    if (i.size() == j.size()) {
      vector<unsigned>::const_iterator i_it = i.begin();
      vector<unsigned>::const_iterator j_it = j.begin();
      while ((i_it != i.end()) && (j_it != j.end()) && (*i_it == *j_it)) {
	++i_it;
	++j_it;
      }
      if ((i_it == i.end()) && (j_it == j.end())) {
	return false;
      }
      return (*i_it < *j_it);
    }
    return (i.size() < j.size());
  }
};

struct compBitsetClusters {
  bool operator() (const boost::dynamic_bitset<unsigned long>& i, const boost::dynamic_bitset<unsigned long>& j) const {
    unsigned i_count = i.count();
    unsigned j_count = j.count();
    if (i_count == j_count) {
      return (i < j);
    }
    return (i_count < j_count);
  }
};

bool compClustersDescending(const vector<unsigned>& i, const vector <unsigned>& j) {
  if (i.size() == j.size()) {
    vector<unsigned>::const_iterator i_it = i.begin();
    vector<unsigned>::const_iterator j_it = j.begin();
    while ((i_it != i.end()) && (j_it != j.end()) && (*i_it == *j_it)) {
      ++i_it;
      ++j_it;
    }
    if ((i_it == i.end()) && (j_it == j.end())) {
      return false;
    }
    return (*i_it < *j_it);
  }
  return (i.size() > j.size());
}

class validClusterBitset {
 public:
  validClusterBitset(const boost::dynamic_bitset<unsigned long> & cluster, unsigned clustID, double thisWeight) {
    elements = cluster;
    ID = clustID;
    weight = thisWeight;
    numElementsHere = cluster.count();
  }

  validClusterBitset(const vector<unsigned> & cluster, unsigned clustID, double thisWeight, unsigned numNodes) {
    elements = boost::dynamic_bitset<unsigned long>(numNodes);
    for (vector<unsigned>::const_iterator it = cluster.begin(); it != cluster.end(); ++it) {
      elements[*it] = 1;
    }
    ID = clustID;
    weight = thisWeight;
    numElementsHere = elements.count();
  }

  bool operator<(const validClusterBitset & b) const {
    return numElementsHere < b.numElements();
  }

  inline unsigned getID() {
    return ID;
  }

  inline void setID(unsigned newID) {
    ID = newID;
    return;
  }

  inline double getWeight() {
    return weight;
  }

  inline const boost::dynamic_bitset<unsigned long> & getElements() {
    return elements;
  }

  inline unsigned isElement(unsigned i) {
    return elements[i];
  }

  inline unsigned numElements() const {
    return numElementsHere;
  }

  inline void addElement(unsigned newElem) {
    elements[newElem] = 1;
    ++numElementsHere;
    return;
  }

  inline vector<unsigned> getElementsVector() {
    vector<unsigned> result;
    result.reserve(numElements());
    for (unsigned i = elements.find_first(); i < elements.size(); i = elements.find_next(i)) {
      result.push_back(i);
    }
    return result;
  }

 private:
  boost::dynamic_bitset<unsigned long> elements;
  unsigned ID;
  double weight;
  unsigned numElementsHere;
};


class ClusterBitset {
 public:
  ClusterBitset(const boost::dynamic_bitset<unsigned long> & cluster, unsigned long & clustID, unsigned long & thisTrieID, double thisClusterWeight = 0) {
    elements = cluster;
    isClusterNew = true;
    isClusterAddedToExplain = false;
    isClusterUnexplainedCounted = false;
    necessary = false;
    checkedFinal = false;
    //isClusterAddedToNecessary = false;
    ID = clustID;
    trieID = thisTrieID;
    numElementsHere = cluster.count();
    clusterWeight = thisClusterWeight;
    valid = false;
    //numExtensions = 0;
    uniquelyExplainedEdges = 0;
    //wouldCombine = false;
  }

  ClusterBitset() {
    isClusterNew = false;
    isClusterAddedToExplain = false;
    isClusterUnexplainedCounted = false;
    necessary = false;
    checkedFinal = false;
    //isClusterAddedToNecessary = false;
    ID = 0;
    numElementsHere = 0;
    trieID = 0;
    clusterWeight = 0;
    valid = false;
    //numExtensions = 0;
    uniquelyExplainedEdges = 0;
    //wouldCombine = false;
  }

  inline double getWeight() {
    return clusterWeight;
  }

  /*inline void setWouldCombine() {
    wouldCombine = true;
  }

  inline void resetWouldCombine() {
    wouldCombine = false;
  }

  inline bool getWouldCombine() {
    return wouldCombine;
    }*/

  /*inline double getThresh(double alpha) {
    return clusterWeight - alpha*Corrector::correction(numExtensions);
    }*/

  inline double getThresh(double alpha, unsigned numUnexplainedEdges) {
    if (numUnexplainedEdges == 0) {
      return 0;
    }
    return clusterWeight - alpha*Corrector::correction(numUnexplainedEdges);
  }

  /*inline unsigned getNumExtensions() {
    return numExtensions;
    }*/

  /*inline void extended() {
    ++numExtensions;
    }*/

  /*inline void setNumExtensions(unsigned val) {
    numExtensions = val;
    }*/

  inline bool isUnexplainedCounted() {
    return isClusterUnexplainedCounted;
  }

  inline void setUnexplainedCounted() {
    isClusterUnexplainedCounted = true;
  }

  inline void setUnexplainedUncounted() {
    isClusterUnexplainedCounted = false;
  }

  inline unsigned getUniquelyExplainedEdges() {
    return uniquelyExplainedEdges; 
  }

  inline void setUniquelyExplainedEdges(unsigned val) {
    uniquelyExplainedEdges = val;
    setUnexplainedCounted();
  }

  inline void setWeight(const double & weight) {
    clusterWeight = weight;
    //historicalWeights.push_back(make_pair(weight, numElementsHere));
  }

  inline vector<pair<double,unsigned> >::iterator historicalWeightsBegin() {
    return historicalWeights.begin();
  }

  inline vector<pair<double,unsigned> >::iterator historicalWeightsEnd() {
    return historicalWeights.end();
  }

  inline void clearHistoricalWeights() {
    historicalWeights.clear();
  }

  inline unsigned long getID() {
    return ID;
  }

  inline unsigned long getTrieID() {
    return trieID;
  }

  inline void setTrieID(unsigned long newTrieID) {
    trieID = newTrieID;
  }

  inline void addElement(unsigned newElem) {
    setNew();
    elements[newElem] = 1;
    ++numElementsHere;
    if (elementsVector.size()) {
      Utils::insertInOrder(elementsVector, newElem);
    }
    return;
  }
  
  inline const boost::dynamic_bitset<unsigned long> & getElements() {
    return elements;
  }

  inline const vector<unsigned> & getElementsVector() {
    if (elementsVector.size()) {
      return elementsVector;
    }
    elementsVector.reserve(numElements());
    for (unsigned i = 0; i < elements.size(); ++i) {
      if (elements[i] == true) {
	elementsVector.push_back(i);
      }
    }
    return elementsVector;
  }

  inline unsigned numElements() const {
    return numElementsHere;
  }

  inline bool isNew() {
    return isClusterNew;
  }

  inline bool isAddedToExplain() {
    return isClusterAddedToExplain;
  }

  inline void setAddedToExplain() {
    isClusterAddedToExplain = true;
  }

  inline void setRemovedFromExplain() {
    isClusterAddedToExplain = false;
  }

  inline void setOld() {
    isClusterNew = false;
  }

  inline void setNew() {
    isClusterNew = true;
    necessary = false;
    checkedFinal = false;
    isClusterUnexplainedCounted= false;
  }

  inline bool isElement(unsigned elemID) {
    return elements.test(elemID);
  }

  inline unsigned size() {
    return elements.size();
  }

  inline void setValid() { 
    valid = true;
    necessary = true;
    checkedFinal = true;
  }

  inline void setInvalid() {
    /*if (valid) {
      numExtensions = 0;
      }*/
    valid = false;
  }

  inline bool isValid() {
    return valid;
    /*if (valid) {
      return 1;
    }
    return 0;*/
  }

  inline bool wasNecessary() {
    return necessary;
  }

  inline void setNecessary() {
    necessary = true;
  }

  inline bool wasCheckedFinal() {
    return checkedFinal;
  }

  inline void setCheckedFinal() {
    checkedFinal = true;
  }

 private:
  boost::dynamic_bitset<unsigned long> elements;
  bool isClusterNew;
  bool isClusterAddedToExplain;
  bool isClusterUnexplainedCounted;
  bool necessary;
  bool checkedFinal;
  unsigned long ID;
  unsigned numElementsHere;
  unsigned long trieID;
  vector<unsigned> elementsVector;
  double clusterWeight;
  vector<pair<double,unsigned> > historicalWeights;
  bool valid;
  //unsigned numExtensions;
  unsigned uniquelyExplainedEdges;
  //bool wouldCombine;
};


struct compListClusterBitsetIterator {
  bool operator() (const list<ClusterBitset>::iterator & i, const list<ClusterBitset>::iterator & j) {
    return (i->getID() < j->getID());
  }
};

bool compClusterBitsetObjs(const ClusterBitset & i, const ClusterBitset & j) {
  return (i.numElements() > j.numElements());
}

class clusterMapClass {
 public:
  clusterMapClass() {};
  
  inline bool addKey(const boost::dynamic_bitset<unsigned long> & key) {
    return Utils::insertInOrder(clustersVec, key);
  }

  inline bool deleteKey(const boost::dynamic_bitset<unsigned long> & key) {
    return Utils::eraseInOrder(clustersVec, key);
  }

  inline bool keyExists(const boost::dynamic_bitset<unsigned long> & key) {
    return Utils::elementExists(clustersVec, key);
  }

 private:
  vector<boost::dynamic_bitset<unsigned long> > clustersVec;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


class currentClusterClassBitset {
 public:
  currentClusterClassBitset(unsigned numNodesHere, unsigned numBlocks, double FirstWeight = 1, double thisAlpha = 0) { // : clusterTrie(Trie(numNodesHere*1000)) {
    //currentClusterClassBitset(unsigned numNodesHere, unsigned numBlocks, double FirstWeight = 1) {
    numNodes = numNodesHere;
    nodesToClusters = vector<vector<unsigned long> >(numNodes, vector<unsigned long>());
    nextID = 0;
    clustersDeleted = 0;
    clustersAdded = 0;
    clustersExtended = 0;
    newClusts = 0;
    firstWeight = FirstWeight;
    curWeight = firstWeight;
    maxNewWeight = firstWeight;
    minWeightAdded = firstWeight;
    combiningNow = false;
    alpha = thisAlpha;
    
    edgesToClusters = vector<vector<unsigned> >(numNodes, vector<unsigned>(numNodes, 0));
    isEdgeExplained = vector<vector<char> >(numNodes, vector<char>(numNodes, 0));
  };

  /*inline void printCluster(const boost::dynamic_bitset<unsigned long> & cluster) {
    for (unsigned i = 0; i < cluster.size(); ++i) {
      if (cluster[i]) {
	cout << i << ",";
      }
    }
    }*/

  inline void setCurWeight(const double & weight) {
    if (weight < minWeightAdded) {
      minWeightAdded = weight;
    }
    curWeight = weight;
  }

  /*inline unsigned numTrieNodes() {
    return clusterTrie.numTrieNodes();
    }*/
  
  inline void setEdgeExplained(const unsigned & i, const unsigned & j) {
    isEdgeExplained[i][j] = 1;
  }

  inline bool isThisEdgeExplained(const unsigned & i, const unsigned & j) {
    return isEdgeExplained[i][j];
  }

  /*inline bool clusterExists(const boost::dynamic_bitset<unsigned long> & cluster) {
    return clusterTrie.keyExists(cluster);
    }*/

  inline bool addClusterToExplainEdge(const unsigned & edgeEnd1, const unsigned & edgeEnd2, const unsigned long & clustID) {
    ++edgesToClusters[edgeEnd1][edgeEnd2];
    return true;
  }

  inline bool addClusterToExplanation(const vector<unsigned> & cluster, const unsigned long & clustID) {
    bool validClust = false;
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	validClust = addClusterToExplainEdge(*it1,*it2,clustID);
      }
    }
    currentClusters[clustID].setAddedToExplain();
    return validClust;
  }

  inline void extendExplanation(unsigned long clustID, unsigned nodeToAdd) {
    unsigned i = currentClusters[clustID].getElements().find_first();
    while (i < nodeToAdd) {
      addClusterToExplainEdge(i, nodeToAdd, clustID);
      i = currentClusters[clustID].getElements().find_next(i);
    }
    while (i < numNodes) {
      addClusterToExplainEdge(nodeToAdd, i, clustID);
      i = currentClusters[clustID].getElements().find_next(i);
    }
    return;
  }

  inline void removeClusterToExplainEdge(const unsigned & edgeEnd1, const unsigned & edgeEnd2, const unsigned long & clustID) {
    --edgesToClusters[edgeEnd1][edgeEnd2];
    return;
  }

  inline void removeClusterFromExplanation(const vector<unsigned> & cluster, const unsigned long & clustID) {
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	removeClusterToExplainEdge(*it1,*it2,clustID);
      }
    }
    currentClusters[clustID].setRemovedFromExplain();
  }

  inline void removeClusterFromExplanation(unsigned long clusterToRemove) {
    removeClusterFromExplanation(currentClusters[clusterToRemove].getElementsVector(), currentClusters[clusterToRemove].getID());
  }

  /*inline*/ unsigned long addCluster(const boost::dynamic_bitset<unsigned long> & newCluster, nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames) {
    unsigned long newClusterTrieID = 0;
    unsigned long newID = 0;
    //if (clusterTrie.addKey(newCluster, newClusterTrieID) == true) {
    if (clusterMap.addKey(newCluster) == true) {
      //unsigned long newID;
      if (openIDs.size() != 0) {
	newID = openIDs.back();
	openIDs.pop_back();
      } else {
	newID = nextID;
	++nextID;
      }
      
      if (currentClusters.size() <= newID) {
	currentClusters.resize(newID+1);
      }
      currentClusters[newID] = ClusterBitset(newCluster, newID, newClusterTrieID, curWeight);
      for (unsigned i = 0; i < newCluster.size() ;  ++i) {
	if (newCluster[i] == true) {
	  Utils::insertInOrder(nodesToClusters[i], newID);
	}
      }
      ++clustersAdded;
      ++newClusts;

      //if (curWeight > minWeightAdded) {
      resetClusterWeight(newID,nodeDistances);
      //if (currentClusters[newID].getWeight() != curWeight) {
      //cout << "RESET WEIGHT CHANGED" << endl;
      //}
      //}
      
      //cout << "# Adding cluster:\t" << newID << "\t" << numElements(newID) << "\t" << curWeight << endl;
      /*cout << "# Adding cluster:\t";
      printCluster(newCluster, nodeIDsToNames);
      cout << "\t" << currentClusters[newID].getWeight() << endl;*/
    }
    return newID;
  }

  inline void printAll(vector<string> & nodeIDsToNames) {
	  /*
    cout << "@ PRINTING CLUSTERS" << endl;
    for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
      if (clustIt->numElements() != 0) {
	cout << "@ ";
	printCluster(clustIt->getElements(), nodeIDsToNames);
	cout << endl;
      }
    }
    cout << "@ DONE PRINTING CLUSTERS" << endl;
    */
  }

  inline void resetAllUnexplained() {
    for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
      if ((clustIt->numElements() != 0) && (!clustIt->isValid())) {
	clustIt->setUnexplainedUncounted();
      }
    }
  }

  inline unsigned setNumUniquelyUnexplainedEdges(unsigned long id) {
    unsigned ret = 0;
    const vector<unsigned> cluster = currentClusters[id].getElementsVector();
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	if (!isThisEdgeExplained(*it1,*it2)) {
	  if (edgesToClusters[*it1][*it2] == 1) {
	    ++ret;
	  }
	}
      }
    }
    currentClusters[id].setUniquelyExplainedEdges(ret);
    return ret;
  }

  inline void setAllNumUniquelyExplained() {
    for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin();
	 clustIt != currentClusters.end(); ++clustIt) {
      if ((clustIt->numElements() != 0) && (!clustIt->isValid())) {
	setNumUniquelyUnexplainedEdges(clustIt->getID());
      }
    }
  }
  
  /*inline*/ void setClusterValid(const boost::dynamic_bitset<unsigned long> & cluster) {
    vector<unsigned> clusterElems;
    clusterElems.reserve(cluster.size());
    for (unsigned i = cluster.find_first(); i < cluster.size(); i = cluster.find_next(i)) {
      clusterElems.push_back(i);
    }
    for (unsigned i = 0; i < clusterElems.size()-1;  ++i) {
      for (unsigned j = i+1; j < clusterElems.size(); ++j) {
	setEdgeExplained(clusterElems[i],clusterElems[j]);
      }
    }
    return;
  }

  inline void setClusterValid(unsigned long clusterID) {
    setNumUniquelyUnexplainedEdges(clusterID);
    setClusterValid(getElements(clusterID));
    currentClusters[clusterID].setValid();
    return;
  }

  inline void setNecessary(unsigned long clusterID) {
    currentClusters[clusterID].setNecessary();
    return;
  }

  inline bool wasNecessary(unsigned long clusterID) {
    return currentClusters[clusterID].wasNecessary();
  }

  inline void setCheckedFinal(unsigned long clusterID) {
    currentClusters[clusterID].setCheckedFinal();
    return;
  }

  inline bool wasCheckedFinal(unsigned long clusterID) {
    return currentClusters[clusterID].wasCheckedFinal();
  }

  inline double getMaxThresh(unsigned long id) {
    //return currentClusters[id].getThresh(alpha);
    return currentClusters[id].getThresh(alpha, getNumUnexplainedEdges(id));
  }

  inline double getThresh(unsigned long id) {
    //return currentClusters[id].getThresh(alpha);
    if (!currentClusters[id].isValid() && !currentClusters[id].isUnexplainedCounted()) {
      setNumUniquelyUnexplainedEdges(id);
    }
    return currentClusters[id].getThresh(alpha, getNumUniquelyUnexplainedEdges(id));
  }

  inline double getMaxThresh() {
    unsigned numFound = 0;
    nextThreshold = 0;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->isNew()) {
	++numFound;
	//double clustThresh = clustIt->getThresh(alpha);
	double clustThresh = getMaxThresh(clustIt->getID());
	if (clustThresh > nextThreshold) {
	  nextThreshold = clustThresh;
	}
      }
      if (numFound == numNew()) {
	return nextThreshold;
      }
    }
    return nextThreshold;
  }

  inline double getNextThresh() {
    unsigned numFound = 0;
    nextThreshold = 0;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->isNew()) {
	++numFound;
	//double clustThresh = clustIt->getThresh(alpha);
	//double clustThresh = getMaxThresh(clustIt->getID());
	double clustThresh = getThresh(clustIt->getID());
	if (clustThresh > nextThreshold) {
	  nextThreshold = clustThresh;
	}
      }
      if (numFound == numNew()) {
	return nextThreshold;
      }
    }
    return nextThreshold;
  }

  inline double getMaxNewWeight() {
    unsigned numFound = 0;
    maxNewWeight = curWeight;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->isNew()) {
	++numFound;
	if (clustIt->getWeight() > maxNewWeight) {
	  maxNewWeight = clustIt->getWeight();
	}
      }
      if (numFound == numNew()) {
	return maxNewWeight;
      }
    }
    return maxNewWeight;
  }

  inline void resetMaxWeight() {
    unsigned numFound = 0;
    maxNewWeight = curWeight;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->isNew()) {
	++numFound;
	if (clustIt->getWeight() > maxNewWeight) {
	  maxNewWeight = clustIt->getWeight();
	}
      }
      if (numFound == numNew()) {
	//cout << "maxNewWeight: " << maxNewWeight << endl;
	return;
      }
    }
    //cout << "WHY? maxNewWeight: " << maxNewWeight << endl;
    return;
  }

  bool isMinNodeDegreeMet(unsigned cluster1, unsigned cluster2, graph_undirected_bitset & clusterGraph, double density) {
    boost::dynamic_bitset<unsigned long> combination = getElements(cluster1);
    combination |= getElements(cluster2);
    unsigned numCombined = combination.count();
    double denom = numCombined - 1;
    unsigned allButOne = numCombined - 2;
    unsigned numChecked = 0;
    //for (unsigned i = 0; i < combination.size(); ++i) {
    //if (combination[i] == 1) {
    for (unsigned i = combination.find_first(); i < combination.size(); i = combination.find_next(i)) {
      boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
      interactorsInCombo &= clusterGraph.getInteractors(i);
      unsigned numInteractorsInCombo = interactorsInCombo.count();
      if ((numInteractorsInCombo / denom) < density) {
	if (numInteractorsInCombo < allButOne) {
	  return false;
	}
      }
      ++numChecked;
      if (numChecked == numCombined) {
	return true;
      }
    }
    //}
    return true;
  }

  inline bool wouldClusterCombine(const unsigned long i, graph_undirected_bitset & clusterGraph, double density) {
    for (unsigned long j = 0; j < maxClusterID(); ++j) {
      if ((j != i) && (numElements(j) != 0)) {
	if (isMinNodeDegreeMet(i,j,clusterGraph,density)) {
	  return true;
	}
      }
    }
    return false;
  }

  
  inline void deleteCluster(const unsigned long & clusterToDelete, vector<string> & nodeIDsToNames, bool printClusterInfo = true) {
    /*cout << "# Deleting cluster:\t";
    printCluster(currentClusters[clusterToDelete].getElements(), nodeIDsToNames);
    cout << endl;*/

    //if (!combiningNow) {
    //if (currentClusters[clusterToDelete].wasNecessary() && !(currentClusters[clusterToDelete].wasCheckedFinal() && !currentClusters[clusterToDelete].isValid())) {
    //if (currentClusters[clusterToDelete].wasNecessary()) {
    //if (printClusterInfo || currentClusters[clusterToDelete].wasNecessary()) {
    

    //if (printClusterInfo && (!combiningNow || (currentClusters[clusterToDelete].wasCheckedFinal() && currentClusters[clusterToDelete].isValid()))) {

    // REMOVED TEMPORARILY
    /*if ((printClusterInfo && /*!combiningNow &&*//* !currentClusters[clusterToDelete].wasCheckedFinal()) || (currentClusters[clusterToDelete].isValid())) {
      cout << "#$$\t" << numElements(clusterToDelete) << "\t" << getClusterWeight(clusterToDelete) - getMinWeightAdded() << "\t" << currentClusters[clusterToDelete].isValid() << "\t" << currentClusters[clusterToDelete].wasNecessary() << "\t" << getClusterWeight(clusterToDelete) << "\t" << getMinWeightAdded() << "\t" << getNumUniquelyUnexplainedEdges(clusterToDelete) << "\t" << getNumUnexplainedEdges(clusterToDelete) << "\t" << combiningNow << "\t" << currentClusters[clusterToDelete].getWouldCombine() << '\t';
      printCluster(currentClusters[clusterToDelete].getElements(), nodeIDsToNames);
      cout << endl;
      }*/

    if (currentClusters[clusterToDelete].isAddedToExplain()) {
      removeClusterFromExplanation(clusterToDelete);
    }
    openIDs.push_back(clusterToDelete);
    vector<ClusterBitset>::iterator clusterToDelete_it = currentClusters.begin();
    clusterToDelete_it += clusterToDelete;
    // it iterates through the nodes in the cluster to be deleted.  This allows us
    // to remove the now deleted cluster from the nodeToClusters structure
    for (unsigned i = 0; i < numNodes; ++i) {
      if (clusterToDelete_it->isElement(i)) {
	Utils::eraseInOrder(nodesToClusters[i], clusterToDelete);
      }
    }
    //clusterTrie.deleteKeyFromNode(clusterToDelete_it->getTrieID());
    clusterMap.deleteKey(clusterToDelete_it->getElements());
    ++clustersDeleted;
    if (clusterToDelete_it->isNew()) {
      --newClusts;
    }
    *clusterToDelete_it = ClusterBitset();

    //cout << "# Deleting cluster:\t" << clusterToDelete << endl;
    return;
  }

  /*inline*/ void deleteClusters(vector<unsigned long> & clustersToDelete, vector<string> & nodeIDsToNames, graph_undirected_bitset & clusterGraph) {

    // PRINTING STATS
    /*cout << "# DELETING MULTIPLE" << endl;*/
    /*for (vector<unsigned long>::iterator clustersToDelete_it = clustersToDelete.begin();
	 clustersToDelete_it != clustersToDelete.end(); ++clustersToDelete_it) {
      if ((//!combiningNow// && !currentClusters[*clustersToDelete_it].wasCheckedFinal()) || (currentClusters[*clustersToDelete_it].isValid())) {
	cout << "#$$\t" << numElements(*clustersToDelete_it) << "\t" << getClusterWeight(*clustersToDelete_it) - getMinWeightAdded() << "\t" << currentClusters[*clustersToDelete_it].isValid() << "\t" << currentClusters[*clustersToDelete_it].wasNecessary() << "\t" << getClusterWeight(*clustersToDelete_it) << "\t" << getMinWeightAdded() << "\t" << getNumUniquelyUnexplainedEdges(*clustersToDelete_it) << "\t" << getNumUnexplainedEdges(*clustersToDelete_it) << "\t" << combiningNow << "\t" <<  currentClusters[*clustersToDelete_it].getWouldCombine() << "\t";
	printCluster(currentClusters[*clustersToDelete_it].getElements(), nodeIDsToNames);
	cout << endl;
      }
      }*/

    for (vector<unsigned long>::iterator clustersToDelete_it = clustersToDelete.begin();
	 clustersToDelete_it != clustersToDelete.end(); ++clustersToDelete_it) {
      deleteCluster(*clustersToDelete_it, nodeIDsToNames);
    }

    return;
  }

  inline unsigned numNew() {
    return newClusts;
  }

  /*inline unsigned getNumExtensions(unsigned long id) {
    return currentClusters[id].getNumExtensions();
    }*/

  /*inline void extended(unsigned long id) {
    currentClusters[id].extended();
    }*/

  /*inline void setNumExtensions(unsigned long id, unsigned val) {
    currentClusters[id].setNumExtensions(val);
    }*/
  
  /*inline*/ void extendClusters(vector<unsigned long> & clustersToExtend, unsigned nodeToAdd, nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames, graph_undirected_bitset & clusterGraph) {
    for (vector<unsigned long>::iterator clustersToExtend_it = clustersToExtend.begin();
	 clustersToExtend_it != clustersToExtend.end(); ++clustersToExtend_it) {
      
      vector<ClusterBitset>::iterator clusterToExtend_it = currentClusters.begin();
      clusterToExtend_it += *clustersToExtend_it;

      /*cout << "# Deleting by extension of cluster:\t";
      printCluster(clusterToExtend_it->getElements(), nodeIDsToNames);
      cout << "\t" << curWeight << endl;*/

      //if (!combiningNow || (clusterToExtend_it->wasCheckedFinal() && clusterToExtend_it->isValid())) {

      //REMOVED TEMPORARILY
      /*if ((/*!combiningNow && *//*!clusterToExtend_it->wasCheckedFinal()) || clusterToExtend_it->isValid()) {
	
	
	//if (clusterToExtend_it->wasNecessary() && !(clusterToExtend_it->wasCheckedFinal() && !clusterToExtend_it->isValid())) {
	//if (clusterToExtend_it->wasNecessary()) {
	cout << "#$$\t" << clusterToExtend_it->numElements() << "\t" << clusterToExtend_it->getWeight() - getMinWeightAdded() << "\t" << clusterToExtend_it->isValid() << "\t" << clusterToExtend_it->wasNecessary() << "\t" << clusterToExtend_it->getWeight() << "\t" << getMinWeightAdded() << "\t" << getNumUniquelyUnexplainedEdges(*clustersToExtend_it) << "\t" << getNumUnexplainedEdges(*clustersToExtend_it) << "\t" << combiningNow << "\t" <<  clusterToExtend_it->getWouldCombine() << '\t';;
	printCluster(clusterToExtend_it->getElements(), nodeIDsToNames);
	cout << endl;
	}*/
      
      clusterToExtend_it->setInvalid();
      //clusterToExtend_it->extended();
      
      if (clusterToExtend_it->isAddedToExplain()) {
	extendExplanation(*clustersToExtend_it, nodeToAdd);
      }

      //cout << "Extending cluster " << (*clustersToExtend_it)->getID() << endl;
      // Erase the cluster to be extended from the clusterSet, extend it, and then put it back in
      //clusterTrie.deleteKeyFromNode(clusterToExtend_it->getTrieID());
      clusterMap.deleteKey(clusterToExtend_it->getElements());
      if (!(clusterToExtend_it->isNew())) {
	++newClusts;
      }
      clusterToExtend_it->addElement(nodeToAdd);
      unsigned long newClusterTrieID = 0;
      //clusterTrie.addKey(clusterToExtend_it->getElements(), newClusterTrieID);
      clusterMap.addKey(clusterToExtend_it->getElements());
      clusterToExtend_it->setTrieID(newClusterTrieID);
      if (curWeight < clusterToExtend_it->getWeight()) {
	clusterToExtend_it->setWeight(curWeight);
      }
      
      //cout << "# Extending cluster:\t" << *clustersToExtend_it << "\t" << clusterToExtend_it->numElements() << "\t" << curWeight << endl;
      /*cout << "# Extending cluster:\t";
      printCluster(clusterToExtend_it->getElements());
      cout << "\t" << curWeight << endl;*/


      // Add the cluster to be extended to the set of clusters containing the nodeToAdd
      Utils::insertInOrder(nodesToClusters[nodeToAdd], *clustersToExtend_it);
      ++clustersExtended;

      /*cout << "# Added by extension:\t";
      printCluster(clusterToExtend_it->getElements(), nodeIDsToNames);
      cout << "\t" << curWeight << endl;*/

      resetClusterWeight(*clustersToExtend_it, nodeDistances);
    }
  }
  
  inline unsigned numClustersWithNode(unsigned nodeID) {
    return nodesToClusters[nodeID].size();
  }
  
  inline vector<ClusterBitset>::iterator clustersBegin() {
    return currentClusters.begin();
  }
  
  inline vector<ClusterBitset>::iterator clustersEnd() {
    return currentClusters.end();
  }
  
  inline unsigned numCurrentClusters() {
    return currentClusters.size() - openIDs.size();
  }

  inline vector<unsigned long>::iterator clustersWithNodeBegin(unsigned nodeID) {
    return nodesToClusters[nodeID].begin();
  }

  inline vector<unsigned long>::iterator clustersWithNodeEnd(unsigned nodeID) {
    return nodesToClusters[nodeID].end();
  }

  /*void sortNewClusters(vector<unsigned long> & sortedNewClusters) {
    vector<pair<unsigned long, pair<unsigned, unsigned> > > newClustersAndCounts;
    newClustersAndCounts.reserve(numCurrentClusters());
    unsigned numFound = 0;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->numElements() != 0) {
	newClustersAndCounts.push_back(make_pair(clustIt->getID(), make_pair(clustIt->numElements(), clustIt->getNumExtensions())));
	++numFound;
      }
    }
    sort(newClustersAndCounts.begin(), newClustersAndCounts.end(), compNewClustersAndCounts);
    
    sortedNewClusters.reserve(numCurrentClusters());
    for (vector<pair<unsigned long, pair<unsigned, unsigned> > >::iterator it = newClustersAndCounts.begin();
	 it != newClustersAndCounts.end(); ++it) {
      sortedNewClusters.push_back(it->first);
    }
    return;
    }*/

  /*inline*/ void sortNewClusters(vector<unsigned long> & sortedNewClusters) {
    vector<pair<unsigned long, unsigned> > newClustersAndCounts;
    newClustersAndCounts.reserve(numCurrentClusters());
    unsigned numFound = 0;
    for (vector<ClusterBitset>::reverse_iterator clustIt = currentClusters.rbegin();
	 clustIt != currentClusters.rend(); ++clustIt) {
      if (clustIt->numElements() != 0) {
	newClustersAndCounts.push_back(make_pair(clustIt->getID(), clustIt->numElements()));
	++numFound;
      }
    }
    sort(newClustersAndCounts.begin(), newClustersAndCounts.end(), compPairSecondAscending);
    
    sortedNewClusters.reserve(numCurrentClusters());
    for (vector<pair<unsigned long, unsigned> >::iterator it = newClustersAndCounts.begin();
	 it != newClustersAndCounts.end(); ++it) {
      sortedNewClusters.push_back(it->first);
    }
    return;
  }

  inline void clearEdgesToClusters() {
    for (unsigned i = 0; i < numNodes; ++i) {
      for (unsigned j = (i+1); j < numNodes; ++j) {
	edgesToClusters[i][j] = 0;
      }
    }
    return;
  }

  inline void addClustersToExplanations(vector<unsigned long> & sortedNewClusters) {
    for (vector<unsigned long>::iterator newClustIt = sortedNewClusters.begin();
	 newClustIt != sortedNewClusters.end(); ++newClustIt) {
      if (!currentClusters[*newClustIt].isAddedToExplain()) {
	addClusterToExplanation(currentClusters[*newClustIt].getElementsVector(), *newClustIt);
      }
    }
    return;
  }
  
  
  /*inline*/ void prepareForValidityCheck(vector<unsigned long> & sortedNewClusters) {
    sortNewClusters(sortedNewClusters);
    addClustersToExplanations(sortedNewClusters);
    return;
  }
  
  bool isLargestExplainer(unsigned edgeEnd1, unsigned edgeEnd2, unsigned long id, vector<char> & idsChecked) {
    unsigned smallest = edgeEnd1;
    unsigned other = edgeEnd2;
    if (numClustersWithNode(edgeEnd2) < numClustersWithNode(edgeEnd1)) {
      smallest = edgeEnd2;
      other = edgeEnd1;
    }
    for (vector<unsigned long>::iterator nodeToClustersIt = clustersWithNodeBegin(smallest);
         nodeToClustersIt != clustersWithNodeEnd(smallest); ++nodeToClustersIt) {
      if (!idsChecked[*nodeToClustersIt] && currentClusters[*nodeToClustersIt].isElement(other)) {
	return false;
      }
    }
    return true;
  }
  
  // Will set isNecessary to true if is only cluster covering an edge
  // Will return true if cluster is only cluster covering an edge which is not previously explained
  inline bool checkClusterFinalValidity(unsigned long id, bool & isNecessary, vector<char> & idsChecked, bool checkForFinal = true) {
    isNecessary = false;
    vector<pair<unsigned, unsigned> > unexplainedEdges;
    const vector<unsigned> cluster = currentClusters[id].getElementsVector();
    unexplainedEdges.reserve(cluster.size()*(cluster.size()-1)/2);
    idsChecked[id] = 1;

    if (!checkForFinal) {
      for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
	vector<unsigned>::const_iterator it2 = it1;
	++it2;
	for ( ; it2 != cluster.end(); ++it2) {
	  if (edgesToClusters[*it1][*it2] == 1) {
	    isNecessary = true;
	    return false;
	  }
	}
      }
    } else {
      for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
	vector<unsigned>::const_iterator it2 = it1;
	++it2;
	for ( ; it2 != cluster.end(); ++it2) {
	  if (edgesToClusters[*it1][*it2] == 1) {
	    isNecessary = true;
	    if (!isThisEdgeExplained(*it1,*it2)) {
	      return true;
	    }
	  }
	  if (!isThisEdgeExplained(*it1,*it2)) {
	    unexplainedEdges.push_back(make_pair(*it1,*it2));
	  }
	}
      }
      if (checkForFinal && isNecessary) {
	//cout << "Necessary but not yet valid" << endl;
	for (vector<pair<unsigned, unsigned> >::iterator unexplainedIt = unexplainedEdges.begin();
	     unexplainedIt != unexplainedEdges.end(); ++unexplainedIt) {
	  if (isLargestExplainer(unexplainedIt->first, unexplainedIt->second, id, idsChecked)) {
	    //cout << "Now valid" << endl;
	    //cout << id << " explains " << unexplainedIt->first << " " << unexplainedIt->second << " via largest possible" << endl;
	    return true;
	  }
	}
	//cout << "Not valid" << endl;
      }
    }
    return false;
  }

  inline unsigned getNumUnexplainedEdges(unsigned long id) {
    unsigned ret = 0;
    const vector<unsigned> cluster = currentClusters[id].getElementsVector();
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	if (!isThisEdgeExplained(*it1,*it2)) {
	  //if (edgesToClusters[*it1][*it2] == 1) {
	    ++ret;
	    //}
	}
      }
    }
    return ret;
  }

  inline unsigned getNumUniquelyUnexplainedEdges(unsigned long id) {
    return currentClusters[id].getUniquelyExplainedEdges();
  }

  inline bool checkClusterValidity(unsigned long id) {
    // 0 is equivalent to false
    const vector<unsigned> cluster = currentClusters[id].getElementsVector();
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	if (edgesToClusters[*it1][*it2] == 1) {
	  setNecessary(id);
	  return true;
	}
      }
    }
    return false;
  }

  inline const boost::dynamic_bitset<unsigned long> & getElements(unsigned long id) {
    return currentClusters[id].getElements();
  }

  inline unsigned numElements(unsigned long id) {
    return currentClusters[id].numElements();
  }

  inline void setOld(unsigned long id) {
    if (currentClusters[id].isNew()){
      --newClusts;
    }
    return currentClusters[id].setOld();
  }

  inline void setNew(unsigned long id) {
    if (!currentClusters[id].isNew()) {
      ++newClusts;
    }
    return currentClusters[id].setNew();
  }

  inline bool isNew(unsigned long id) {
    return currentClusters[id].isNew();
  }

  inline double getClusterWeight(unsigned long id) {
    return currentClusters[id].getWeight();
  }

  inline unsigned long maxClusterID() {
    return currentClusters.size();
  }

  inline double getCurWeight() {
    return curWeight;
  }

  inline double getMinWeightAdded() {
    return minWeightAdded;
  }

  void resetClusterWeight(unsigned long id, nodeDistanceObject & nodeDistances) {
    double minWeight = firstWeight;
    const vector<unsigned> cluster = currentClusters[id].getElementsVector();
    for (vector<unsigned>::const_iterator it1 = cluster.begin(); it1 != cluster.end(); ++it1) {
      vector<unsigned>::const_iterator it2 = it1;
      ++it2;
      for ( ; it2 != cluster.end(); ++it2) {
	//if (!isThisEdgeExplained(*it1,*it2) && (edgesToClusters[*it1][*it2] == 1)) {
	double thisDistance = nodeDistances.getDistance(*it1,*it2);
	if (thisDistance < minWeight) {
	  minWeight = thisDistance;
	}
      }
    }
    currentClusters[id].setWeight(minWeight);
    return;
  }

  /*inline void setWouldCombine(vector<unsigned long> clusts, graph_undirected_bitset & clusterGraph, double density) {
    for (vector<unsigned long>::iterator it = clusts.begin(); it != clusts.end(); ++it) {
      if (wouldClusterCombine(*it, clusterGraph, density)) {
	currentClusters[*it].setWouldCombine();
      }
    }
    }*/

  /*inline void resetWouldCombine(vector<unsigned long> clusts) {
    for (vector<unsigned long>::iterator it = clusts.begin(); it != clusts.end(); ++it) {
      currentClusters[*it].resetWouldCombine();
    }
    } */   

  bool isClusterMaximal(unsigned long id, graph_undirected_bitset & clusterGraph, nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames) {
    bool isValid = checkClusterValidity(id);
    //bool extended = false;
    unsigned nextElem = currentClusters[id].getElements().find_first();
    for (unsigned i = 0; i < numNodes; ++i) {
      if (i == nextElem) {
	nextElem = currentClusters[id].getElements().find_next(nextElem);
      } else if (currentClusters[id].getElements().is_subset_of(clusterGraph.getInteractors(i))) {
	if (isValid) {
	  //cout << "# Extending necessary non-maximal cluster" << endl;
	  vector<unsigned long> clusterToExtend(1,id);
	  //cout << "# ";
	  //printCluster(currentClusters[id].getElements(), nodeIDsToNames);
	  //cout << "\t";
	  //setWouldCombine(clusterToExtend, clusterGraph);
	  extendClusters(clusterToExtend, i, nodeDistances, nodeIDsToNames, clusterGraph);
	  //cout << nodeIDsToNames[i] << "\t";
	  //printCluster(currentClusters[id].getElements(), nodeIDsToNames);
	  //cout << endl;
	  //currentClusters[clusterToExtend[0]].resetWouldCombine();
	  //extended = true;
	} /*else {
	  return false;
	  }*/
      }
    }
    /*if (extended) {
      resetClusterWeight(id, nodeDistances);
      }*/
    return true;
  }


  void removeNonMaximalClusters(graph_undirected_bitset & clusterGraph, nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames) {
    //cout << "# Removing non-maximal" << endl;
    //bool alreadyCombining = combiningNow;
    //combining();
    vector<unsigned long> newClustersSorted;
    prepareForValidityCheck(newClustersSorted);
    for (vector<unsigned long>::iterator newClusterIt = newClustersSorted.begin();
	 newClusterIt != newClustersSorted.end(); ++newClusterIt) {
      if (currentClusters[*newClusterIt].isNew()) {
	if (!isClusterMaximal(*newClusterIt, clusterGraph, nodeDistances,nodeIDsToNames)) {
	  deleteCluster(*newClusterIt, nodeIDsToNames, false);
	}
      }
    }
    //if (!alreadyCombining) {
    //doneCombining();
    //}
    //cout << "# Done removing non-maximal" << endl;
  }

  bool checkAllForMaximality(graph_undirected_bitset & clusterGraph, vector<string> & nodeIDsToNames) {
    bool ret = false;
    //cout << "CHECKING FOR MAXIMALITY" << endl;
    //cout << "Num current clusters: " << numCurrentClusters() << endl;
    for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
      if (clustIt->numElements() != 0) {
	bool isMaximal = true;
	unsigned nextElem = clustIt->getElements().find_first();
	for (unsigned i = 0; i < clustIt->size(); ++i) {
	  if (i == nextElem) {
	    nextElem = clustIt->getElements().find_next(nextElem);
	  } else if (clustIt->getElements().is_subset_of(clusterGraph.getInteractors(i))) {
	    isMaximal = false;
	  //  cout << "CLUSTER NOT MAXIMAL" << endl;
	    ret = true;
	    for (unsigned j = clustIt->getElements().find_first(); j < clustIt->size(); j = clustIt->getElements().find_next(j)) {
	      //cout << nodeIDsToNames[j] << " ";
	    }
	    //cout << endl;
	    //cout << "SHOULD ALSO CONTAIN: " << nodeIDsToNames[i] << endl;
	  }
	}
	if (isMaximal) {
	  for (unsigned j = clustIt->getElements().find_first(); j < clustIt->size(); j = clustIt->getElements().find_next(j)) {
	    //cout << nodeIDsToNames[j] << " ";
	  }
	  //cout << endl;
	}
      }
    }
    return ret;
  }

  unsigned numEdgesCovered() {
    graph_undirected_bitset coveredEdges(numNodes);
    for (vector<ClusterBitset>::iterator clustIt = currentClusters.begin(); clustIt != currentClusters.end(); ++clustIt) {
      if (clustIt->numElements() != 0) {
	vector<unsigned> clusterElems = clustIt->getElementsVector();
	for (unsigned i = 0; i < clusterElems.size()-1; ++i) {
	  for (unsigned j = i+1; j < clusterElems.size(); ++j) {
	    if (!coveredEdges.isEdge(clusterElems[i],clusterElems[j])) {
	      coveredEdges.addEdge(clusterElems[i],clusterElems[j]);
	    }
	  }
	}
      }
    }
    return coveredEdges.numEdges();
  }
  
  unsigned clustersAdded;
  unsigned clustersDeleted;
  unsigned clustersExtended;

  inline vector<pair<double,unsigned> >::iterator historicalWeightsBegin(unsigned long id) {
    return currentClusters[id].historicalWeightsBegin();
  }
  
  inline vector<pair<double,unsigned> >::iterator historicalWeightsEnd(unsigned long id) {
    return currentClusters[id].historicalWeightsEnd();
  }

  inline void clearHistoricalWeights(unsigned long id) {
    currentClusters[id].clearHistoricalWeights();
  }

  inline void combining() {
    combiningNow = true;
  }

  inline void doneCombining() {
    combiningNow = false;
  }
  
  inline bool isCombiningNow() {
    return combiningNow;
  }

  inline double getFirstWeight() {
    return firstWeight;
  }

 private:
  //Trie clusterTrie;
  clusterMapClass clusterMap;
  vector<ClusterBitset> currentClusters;
  vector<unsigned long> openIDs;

  // Each node has an ID.  The nodesToClusters object can be accessed by that ID,
  // giving a set of iterators to cluster objects in the currentClusters list
  // which can then be derefenced to get a vector<unsigned> for the cluster
  //vector<set<list<ClusterBitset>::iterator, compListClusterBitsetIterator > > nodesToClusters;
  vector<vector<unsigned long> > nodesToClusters;
  vector<vector<unsigned> > edgesToClusters;
  vector<vector<char> > isEdgeExplained;

  unsigned long numClusters;
  unsigned long nextID;
  unsigned newClusts;
  unsigned numNodes;
  double curWeight;
  double minWeightAdded;
  double firstWeight;
  bool combiningNow;
  double alpha;

  // Maximum weight of any new cluster curently in this list
  double maxNewWeight;
  double nextThreshold;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

namespace dagConstruct {

  inline void setIntersection(const vector<unsigned>& i, const vector <unsigned>& j, vector<unsigned> & intersection) {
    intersection.clear();
    intersection.reserve(i.size());
    intersection.reserve(j.size());
    vector<unsigned>::const_iterator i_it = i.begin();
    vector<unsigned>::const_iterator j_it = j.begin();
    while ((i_it != i.end()) && (j_it != j.end())) {
      if (*i_it == *j_it) {
	intersection.push_back(*i_it);
	++i_it;
	++j_it;
      } else if (*i_it < *j_it) {
	++i_it;
      } else if (*j_it < *i_it) {
	++j_it;
      }
    }
  }

  inline void getNewCluster(const vector<unsigned>& currentCluster, const vector<unsigned>& neighborsOfNewInteractor, 
			    unsigned newInteractor, vector<unsigned> & newCluster) {
    newCluster.clear();
    newCluster.reserve(currentCluster.size() + 1);
    vector<unsigned>::const_iterator currentCluster_it = currentCluster.begin();
    vector<unsigned>::const_iterator neighbors_it = neighborsOfNewInteractor.begin();
    bool addedNewInteractor = false;
    while ((currentCluster_it != currentCluster.end()) && (neighbors_it != neighborsOfNewInteractor.end())) {
      if (*currentCluster_it == *neighbors_it) {
	if (!addedNewInteractor && (*currentCluster_it > newInteractor)) {
	  newCluster.push_back(newInteractor);
	  addedNewInteractor = true;
	}
	newCluster.push_back(*currentCluster_it);
	++currentCluster_it;
	++neighbors_it;
      } else if (*currentCluster_it < *neighbors_it) {
	++currentCluster_it;
      } else if (*neighbors_it < *currentCluster_it) {
	++neighbors_it;
      }
    }
    if (!addedNewInteractor) {
      newCluster.push_back(newInteractor);
    }
    
    return;
  }


  // newCluster = intersection of currentCluster and neighborsOfNewInteractor + newInteractor
  inline void getNewCluster(const boost::dynamic_bitset<unsigned long>& currentCluster, const boost::dynamic_bitset<unsigned long>& neighborsOfNewInteractor, 
			    unsigned newInteractor, boost::dynamic_bitset<unsigned long> & newCluster) {
    newCluster = currentCluster;
    newCluster &= neighborsOfNewInteractor;
    newCluster[newInteractor] = 1;
    return;
  }
  
  bool isClusterAncestor(const vector<unsigned> & cluster, const vector<unsigned> & possibleAncestorCluster, vector<bool> & unaccountedFor) {
    // Keep track of which genes in the ancestor cluster are "accounted for" (i.e. contained in) decsendent
    vector<unsigned> accountedFor;
    accountedFor.reserve(unaccountedFor.size());

    vector<unsigned>::const_iterator clusterIt = cluster.begin();
    vector<unsigned>::const_iterator ancestorIt = possibleAncestorCluster.begin();
    unsigned ancestorPos = 0;
    while ((clusterIt != cluster.end()) && (ancestorIt != possibleAncestorCluster.end())) {
      if (*ancestorIt < *clusterIt) {
	++ancestorIt;
	++ancestorPos;
      } else if (*ancestorIt > *clusterIt) {
	return false;
      } else if (*ancestorIt == *clusterIt) {
	accountedFor.push_back(ancestorPos);
	++ancestorPos;
	++ancestorIt;
	++clusterIt;
      }
    }
    if (clusterIt == cluster.end()) {
      for (vector<unsigned>::iterator accountedForIt = accountedFor.begin(); accountedForIt != accountedFor.end(); ++accountedForIt) {
	unaccountedFor[*accountedForIt] = false;
      }
      return true;
    }
    return false;
  }

  bool isClusterAncestor(const boost::dynamic_bitset<unsigned long> & cluster, const boost::dynamic_bitset<unsigned long> & possibleAncestorCluster, boost::dynamic_bitset<unsigned long> & unaccountedFor) {
    // Keep track of which genes in the ancestor cluster are "accounted for" (i.e. contained in) decsendent
    if (cluster.is_subset_of(possibleAncestorCluster)) {
      unaccountedFor -= cluster;
      return true;
    }
    return false;
  }

  bool isClusterAncestor(const vector<unsigned> & cluster, const vector<unsigned> & possibleAncestorCluster) {
    vector<unsigned>::const_iterator clusterIt = cluster.begin();
    vector<unsigned>::const_iterator ancestorIt = possibleAncestorCluster.begin();
    while ((clusterIt != cluster.end()) && (ancestorIt != possibleAncestorCluster.end())) {
      if (*ancestorIt < *clusterIt) {
	++ancestorIt;
      } else if (*ancestorIt > *clusterIt) {
	return false;
      } else if (*ancestorIt == *clusterIt) {
	++ancestorIt;
	++clusterIt;
      }
    }
    if (clusterIt == cluster.end()) {
      return true;
    }
    return false;
  }

  inline bool checkClusterValidity(const vector<unsigned> & cluster, nodeDistanceObject & nodeDistances, vector<vector<bool> > & pairLCAFound) {
    //cout << "Checking size " << cluster.size() << endl;
    /*cout << "Checking ";
    for (unsigned i = 0; i < cluster.size(); ++i) {
      cout << cluster[i] << " ";
    }
    cout << endl;*/
    bool firstPairFound = false;
    double firstPairDistance;
    for (unsigned i = 0; i < (cluster.size() - 1); ++i) {
      for (unsigned j = (i+1); j < cluster.size(); ++j) {
	if (!pairLCAFound[cluster[i]][cluster[j]]) {
	  //cout << "Pair " << cluster[i] << " " << cluster[j] << " " << nodeDistances.getDistance(cluster[i],cluster[j]) << endl;
	  if (firstPairFound) {
	    if (firstPairDistance != nodeDistances.getDistance(cluster[i],cluster[j])) {
	      return false;
	    }
	  } else {
	    firstPairFound = true;
	    firstPairDistance = nodeDistances.getDistance(cluster[i],cluster[j]);
	  }
	}
      }
    }
    if (firstPairFound) {
      for (unsigned i = 0; i < (cluster.size() - 1); ++i) {
	for (unsigned j = (i+1); j < cluster.size(); ++j) {
	  pairLCAFound[cluster[i]][cluster[j]] = true;
	}
      }
      //clusterWeight = (1 - firstPairDistance) / 2.0;
      return true;
    } else {
      return false;
    }
  }

  inline bool checkClusterValidity(const boost::dynamic_bitset<unsigned long> & cluster, nodeDistanceObject & nodeDistances, vector<vector<bool> > & pairLCAFound) {
    bool firstPairFound = false;
    double firstPairDistance;
    for (unsigned i = 0; i < (cluster.size() - 1); ++i) {
      if (cluster[i] == true) {
	for (unsigned j = (i+1); j < cluster.size(); ++j) {
	  if (cluster[j] == true) {
	    if (!pairLCAFound[i][j]) {
	      //cout << "Pair " << cluster[i] << " " << cluster[j] << " " << nodeDistances.getDistance(cluster[i],cluster[j]) << endl;
	      if (firstPairFound) {
		if (firstPairDistance != nodeDistances.getDistance(i,j)) {
		  return false;
		}
	      } else {
		firstPairFound = true;
		firstPairDistance = nodeDistances.getDistance(i,j);
	      }
	    }
	  }
	}
      }
    }
    if (firstPairFound) {
      for (unsigned i = 0; i < (cluster.size() - 1); ++i) {
	if (cluster[i] == true) {
	  for (unsigned j = (i+1); j < cluster.size(); ++j) {
	    if (cluster[j] == true) {
	      pairLCAFound[i][j] = true;
	    }
	  }
	}
      }
      //clusterWeight = (1 - firstPairDistance) / 2.0;
      return true;
    } else {
      return false;
    }
  }

  inline bool isClusterMaximal(vector<boost::dynamic_bitset<unsigned long> >::reverse_iterator & clustersToAddIt, 
			       vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd, 
			       vector<unsigned long> & clustersToExtend,
			       currentClusterClassBitset & currentClusters) {
    for (vector<unsigned long>::iterator clustersExtended_it = clustersToExtend.begin();
	 clustersExtended_it != clustersToExtend.end(); ++clustersExtended_it) {
      //if ((currentClusters.numElements(*clustersExtended_it) != 0) && 
      if (clustersToAddIt->is_subset_of(currentClusters.getElements(*clustersExtended_it))) {
	return false;
      }
    }
    
    if (clustersToAddIt != newClustersToAdd.rbegin()) {
      vector<boost::dynamic_bitset<unsigned long> >::reverse_iterator possibleMaximalClustersIt = clustersToAddIt;
      --possibleMaximalClustersIt;
      for ( ; possibleMaximalClustersIt != newClustersToAdd.rbegin(); --possibleMaximalClustersIt) {
	if (clustersToAddIt->is_subset_of(*possibleMaximalClustersIt)) {
	  return false;
	}
      }
      if (clustersToAddIt->is_subset_of(*possibleMaximalClustersIt)) {
	return false;
      }
    }
    return true;
  }
  
  // RETURNS WHETHER ANYTHING WAS DELETED
  // (No longer does as follows) RETURNS THE NUMBER OF EXTENSIONS OF THE LARGEST DELETED CLUSTER.
  // THIS INFO IS PASSED TO THE NEW CLUSTER BEING EXTENDED
  /*inline*/ unsigned checkForDelete(currentClusterClassBitset & currentClusters, 
				     unsigned otherNode, unsigned nodeToUse, 
				     boost::dynamic_bitset<unsigned long> & newCluster, 
				     vector<unsigned long> & clustersToDelete) {
    //newCluster[nodeToUse] = 0;
    //if (currentClusters.clusterExists(newCluster)) {
    bool retVal = false;
    //unsigned maxExtendedDelete = 0;
    for (vector<unsigned long>::iterator nodeToClustersIt = currentClusters.clustersWithNodeBegin(otherNode);
	 nodeToClustersIt != currentClusters.clustersWithNodeEnd(otherNode); ++nodeToClustersIt) {
      if (currentClusters.getElements(*nodeToClustersIt).is_proper_subset_of(newCluster)) {
	Utils::insertInOrder(clustersToDelete,*nodeToClustersIt);
	/*if (currentClusters.getNumExtensions(*nodeToClustersIt) > maxExtendedDelete) {
	  maxExtendedDelete = currentClusters.getNumExtensions(*nodeToClustersIt);
	  }*/
	retVal = true;
      }
    }
    return retVal;
    //return maxExtendedDelete;
  }

  // After adding a new edge between nodeToUse and otherNode, use this to get clusters which need to be updated
  inline void getClusterUpdates(unsigned nodeToUse, unsigned otherNode, graph_undirected_bitset & clusterGraph,
				currentClusterClassBitset & currentClusters, set<boost::dynamic_bitset<unsigned long> > & newClustersToAddSet,
				vector<unsigned long> & clustersToDelete,
				vector<unsigned long> & clustersToExtend) {
    unsigned numNodes = clusterGraph.numNodes();
    for (vector<unsigned long>::iterator node1ClustersIt = currentClusters.clustersWithNodeBegin(nodeToUse);
	 node1ClustersIt != currentClusters.clustersWithNodeEnd(nodeToUse); ++node1ClustersIt) {
      if (!currentClusters.getElements(*node1ClustersIt).test(otherNode)) {
	boost::dynamic_bitset<unsigned long> newCluster(numNodes);
	getNewCluster(currentClusters.getElements(*node1ClustersIt), clusterGraph.getInteractors(otherNode), otherNode, newCluster);
	
	unsigned newClusterCount = newCluster.count();
	if (newClusterCount > currentClusters.numElements(*node1ClustersIt)) {
	  clustersToExtend.push_back(*node1ClustersIt);
	  /*unsigned maxExtendedDelete = */ checkForDelete(currentClusters, otherNode, nodeToUse, newCluster, clustersToDelete);
	  /*if (currentClusters.getNumExtensions(*node1ClustersIt) < maxExtendedDelete) {
	  currentClusters.setNumExtensions(*node1ClustersIt, maxExtendedDelete);
	  }*/
	} else if (newClusterCount > 2) {
	  if (newClustersToAddSet.insert(newCluster).second) {
	    //checkForDelete(currentClusters, otherNode, nodeToUse, newCluster, clustersToDelete);
	  }
	}
      }
    }
    return;
  }

  /*inline void printCluster(const boost::dynamic_bitset<unsigned long> & cluster, vector<string> & nodeIDsToNames) {
    for (unsigned i = 0; i < cluster.size(); ++i) {
      if (cluster[i]) {
	cout << nodeIDsToNames[i] << ",";
      }
    }
    }*/

  inline void performValidityCheck(currentClusterClassBitset & currentClusters, graph_undirected_bitset & clusterGraph,
				   nodeDistanceObject & nodeDistances, vector<string> & nodeIDsToNames) {
    //currentClusters.removeNonMaximalClusters(clusterGraph, nodeDistances,nodeIDsToNames);
    vector<unsigned long> newClustersSorted;
    currentClusters.prepareForValidityCheck(newClustersSorted);

    //bool alreadyCombining = currentClusters.isCombiningNow();
    //currentClusters.combining();

    for (vector<unsigned long>::iterator newClusterIt = newClustersSorted.begin();
	 newClusterIt != newClustersSorted.end(); ++newClusterIt) {
      if (!currentClusters.checkClusterValidity(*newClusterIt)) {
	currentClusters.deleteCluster(*newClusterIt, nodeIDsToNames,false);
      }
    }

    currentClusters.resetAllUnexplained();

    //if (!alreadyCombining) {
    //currentClusters.doneCombining();
    //}
    return;
  }

  // Takes an empty vector and returns vector sorted by size (large to small)
  inline void newClustersToAddSetToVector(set<boost::dynamic_bitset<unsigned long> > & newClustersToAddSet, vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd) {
    newClustersToAdd.reserve(newClustersToAddSet.size());
    if (newClustersToAddSet.size() == 1) {
      newClustersToAdd.push_back(*(newClustersToAddSet.begin()));
    } else {
      vector<pair<boost::dynamic_bitset<unsigned long>, unsigned> > newClustersToAddSortingVec;
      newClustersToAddSortingVec.reserve(newClustersToAddSet.size());
      for (set<boost::dynamic_bitset<unsigned long>>::iterator setIt = newClustersToAddSet.begin();
	   setIt != newClustersToAddSet.end(); ++setIt) {
	newClustersToAddSortingVec.push_back(make_pair(*setIt, setIt->count()));
      }
      sort(newClustersToAddSortingVec.begin(), newClustersToAddSortingVec.end(), newClustersComp);
      
      for(vector<pair<boost::dynamic_bitset<unsigned long>, unsigned> >::iterator it = newClustersToAddSortingVec.begin();
	  it != newClustersToAddSortingVec.end(); ++it) {
	newClustersToAdd.push_back(it->first);
      }
    }
    return;
  }
  
  void addEdgeAndUpdateClusters(unsigned firstNode, unsigned secondNode, graph_undirected_bitset & clusterGraph, 
				currentClusterClassBitset & currentClusters, vector<string> & nodeIDsToNames, nodeDistanceObject & nodeDistances) {

    // WHY IS THIS HERE???? CHECK THIS SOON
    //currentClusters.setAllNumUniquelyExplained();
    
    //unsigned firstNode = distanceIt->first.first;
    //unsigned secondNode = distanceIt->first.second;
    //cout << "Adding edge: " << nodeIDsToNames[firstNode] << "\t" << nodeIDsToNames[secondNode] << endl;
    //currentClusters.printAll(nodeIDsToNames);
    
    
    /*unsigned smaller = firstNode;
    unsigned larger = secondNode;
    if (secondNode < firstNode) {
      smaller = secondNode;
      larger = firstNode;
    }
    cout << "Adding edge: " << smaller << "\t" << larger << "\t" << nodeDistances.getDistance(smaller,larger) << endl;*/

    clusterGraph.addEdge(firstNode, secondNode);
    //clusterGraph.printAll(nodeIDsToNames);

    if (clusterGraph.getInteractors(firstNode).intersects(clusterGraph.getInteractors(secondNode))) {
      set<boost::dynamic_bitset<unsigned long> > newClustersToAddSet;
      vector<unsigned long> clustersToDelete;
      vector<unsigned long> clustersToExtend;
      vector<unsigned long> clustersToExtendFirst;
      
      clustersToExtendFirst.reserve(currentClusters.numClustersWithNode(firstNode));
      clustersToExtend.reserve(currentClusters.numClustersWithNode(firstNode) + currentClusters.numClustersWithNode(secondNode));
      
      getClusterUpdates(firstNode, secondNode, clusterGraph, currentClusters, newClustersToAddSet, clustersToDelete, clustersToExtendFirst);
      //currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames);
      //currentClusters.extendClusters(clustersToExtendFirst, secondNode, nodeDistances, nodeIDsToNames);

      //currentClusters.setWouldCombine(clustersToDelete, clusterGraph, density);
      //currentClusters.setWouldCombine(clustersToExtendFirst, clusterGraph, density);
      currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames, clusterGraph);
      currentClusters.extendClusters(clustersToExtendFirst, secondNode, nodeDistances, nodeIDsToNames, clusterGraph);
      //currentClusters.resetWouldCombine(clustersToExtendFirst);

      clustersToDelete.clear();
      getClusterUpdates(secondNode, firstNode, clusterGraph, currentClusters, newClustersToAddSet, clustersToDelete, clustersToExtend);
      //currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames);
      //currentClusters.extendClusters(clustersToExtend, firstNode, nodeDistances, nodeIDsToNames);

      //currentClusters.setWouldCombine(clustersToDelete, clusterGraph, density);
      //currentClusters.setWouldCombine(clustersToExtend, clusterGraph, density);
      currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames, clusterGraph);
      currentClusters.extendClusters(clustersToExtend, firstNode, nodeDistances, nodeIDsToNames, clusterGraph);
      //currentClusters.resetWouldCombine(clustersToExtend);
      
      sort(clustersToDelete.begin(), clustersToDelete.end());
      sort(clustersToExtendFirst.begin(), clustersToExtendFirst.end());
      vector<unsigned long> stillExtended(clustersToExtendFirst.size());
      vector<unsigned long>::iterator newEnd;
      newEnd = set_difference(clustersToExtendFirst.begin(), clustersToExtendFirst.end(), clustersToDelete.begin(), clustersToDelete.end(), stillExtended.begin());
      stillExtended.resize(newEnd - stillExtended.begin());
      clustersToExtend.insert(clustersToExtend.end(), stillExtended.begin(), stillExtended.end());
      clustersToDelete.clear();

      if (newClustersToAddSet.size() > 0) {
	vector<boost::dynamic_bitset<unsigned long> > newClustersToAdd;
	newClustersToAddSetToVector(newClustersToAddSet, newClustersToAdd);
	
	//cout << "Num new clusters: " << newClustersToAdd.size() << endl;
	for (vector<boost::dynamic_bitset<unsigned long>>::reverse_iterator clustersToAddIt = newClustersToAdd.rbegin();
	     clustersToAddIt != newClustersToAdd.rend(); ++clustersToAddIt) {	  
	  if (isClusterMaximal(clustersToAddIt, newClustersToAdd, clustersToExtend, currentClusters)) {
	    /*unsigned maxExtendedDelete = */ checkForDelete(currentClusters, firstNode, secondNode, *clustersToAddIt, clustersToDelete);
	    /*currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames);
	      clustersToDelete.clear();*/
	    /*unsigned newID = */ currentClusters.addCluster(*clustersToAddIt,nodeDistances, nodeIDsToNames);
	    //currentClusters.setNumExtensions(newID, maxExtendedDelete + 1);
	  }
	}
	//currentClusters.setWouldCombine(clustersToDelete, clusterGraph, density);
	currentClusters.deleteClusters(clustersToDelete, nodeIDsToNames, clusterGraph);
      }
    } else {
      boost::dynamic_bitset<unsigned long> newCluster(clusterGraph.numNodes());
      newCluster[firstNode] = 1;
      newCluster[secondNode] = 1;
      currentClusters.addCluster(newCluster,nodeDistances, nodeIDsToNames);
    }
    //currentClusters.printAll(nodeIDsToNames);
  }

  bool isMinNodeDegreeMet(unsigned cluster1, unsigned cluster2, currentClusterClassBitset & currentClusters, graph_undirected_bitset & clusterGraph, double density) {
    boost::dynamic_bitset<unsigned long> combination = currentClusters.getElements(cluster1);
    combination |= currentClusters.getElements(cluster2);
    unsigned numCombined = combination.count();
    double denom = numCombined - 1;
    unsigned allButOne = numCombined - 2;
    unsigned numChecked = 0;
    //for (unsigned i = 0; i < combination.size(); ++i) {
    //if (combination[i] == 1) {
    for (unsigned i = combination.find_first(); i < combination.size(); i = combination.find_next(i)) {
      boost::dynamic_bitset<unsigned long> interactorsInCombo = combination;
      interactorsInCombo &= clusterGraph.getInteractors(i);
      unsigned numInteractorsInCombo = interactorsInCombo.count();
      if ((numInteractorsInCombo / denom) < density) {
	if (numInteractorsInCombo < allButOne) {
	  return false;
	}
      }
      ++numChecked;
      if (numChecked == numCombined) {
	return true;
      }
    }
    //}
    return true;
  }

  // Checks the starting cluster for possible extensions based on the new clusterGraph
  void checkForExtensions(unsigned i, unsigned nextElem, boost::dynamic_bitset<unsigned long> startingCluster, graph_undirected_bitset & clusterGraph, 
			  vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd) {
    for ( ; i < clusterGraph.numNodes(); ++i) {
      if (i == nextElem) {
	nextElem = startingCluster.find_next(nextElem);
      } else if (startingCluster.is_subset_of(clusterGraph.getInteractors(i))) {
	
	// N(i) contains the startingCluster
	boost::dynamic_bitset<unsigned long> extendedCluster = startingCluster;
	extendedCluster[i] = 1;

	// Is there some other element j where j < i and N(j) includes this extended cluster? (Maximality test)
	bool extend = true;
	for (unsigned j = 0; j < i; ++j) {
	  if ((extendedCluster[j] == 0) && (extendedCluster.is_subset_of(clusterGraph.getInteractors(j)))) {

	    // There is some j where j < i and N(j) contains this extended Cluster
	    extend = false;
	    break;
	  }
	}

	if (extend) {
	  checkForExtensions(i+1, nextElem, extendedCluster, clusterGraph, newClustersToAdd);
	}

      }
    }

    // We have reached i == clusterGraph.numNodes().  Add startingCluster to newClustersToAdd
    Utils::insertInOrder(newClustersToAdd, startingCluster);
  }

  // Checks the starting cluster for possible extensions based on the new clusterGraph
  bool checkForExtensions(boost::dynamic_bitset<unsigned long> startingCluster, graph_undirected_bitset & clusterGraph, 
			  vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd) {

    unsigned nextElem = startingCluster.find_first();
    bool foundOne = false;

    for (unsigned i = 0; i < clusterGraph.numNodes(); ++i) {
      if (i == nextElem) {
	nextElem = startingCluster.find_next(nextElem);
      } else if (startingCluster.is_subset_of(clusterGraph.getInteractors(i))) {
	//cout << "# Found extension" << endl;
	// N(i) contains the startingCluster
	boost::dynamic_bitset<unsigned long> extendedCluster = startingCluster;
	extendedCluster[i] = 1;
	
	// Is there some other element j where j < i and N(j) includes this extended cluster? (Maximality test)
	if (foundOne == true) {
	  bool extend = true;
	  for (unsigned j = 0; j < i; ++j) {
	    if ((extendedCluster[j] == 0) && (extendedCluster.is_subset_of(clusterGraph.getInteractors(j)))) {
	      
	      // There is some j where j < i and N(j) contains this extended Cluster
	      extend = false;
	      break;
	    }
	  }
	  
	  if (extend) {
	    checkForExtensions(i+1, nextElem, extendedCluster, clusterGraph, newClustersToAdd);
	  }
	} else if (foundOne == false) {
	  checkForExtensions(i+1, nextElem, extendedCluster, clusterGraph, newClustersToAdd);
	}
	
	foundOne = true;
      }
    }

    return foundOne;
  }

  void updateCliques(graph_undirected_bitset & clusterGraph,
		     vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd,
		     boost::dynamic_bitset<unsigned long> startClust, vector<unsigned> & neighborsOfBoth,
		     vector<char> & neighborsOfBothBits, unsigned i, vector<string> & nodeIDsToNames) {
    if (i == neighborsOfBoth.size()) {
      if (!Utils::insertInOrder(newClustersToAdd, startClust)) {
	//cout << "# Repeated" << endl;
      } /*else {
	cout << "# Adding ";
	printCluster(startClust, nodeIDsToNames);
	cout << endl;
	}*/
      return;
    }
    //printCluster(startClust, nodeIDsToNames);
    //cout << endl;
    if (startClust.is_subset_of(clusterGraph.getInteractors(neighborsOfBoth[i]))) {
      startClust[neighborsOfBoth[i]] = 1;
      updateCliques(clusterGraph, newClustersToAdd, startClust, neighborsOfBoth, neighborsOfBothBits, i + 1, nodeIDsToNames);
    } else {
      updateCliques(clusterGraph, newClustersToAdd, startClust, neighborsOfBoth, neighborsOfBothBits, i + 1, nodeIDsToNames);
      
      // Calculate T[y] = |N(y) intersection with C intersection with N(i)| for y not in C, not i
      vector<unsigned> T(clusterGraph.numNodes(), 0);
      
      boost::dynamic_bitset<unsigned long> nextClust = startClust;
      nextClust &= clusterGraph.getInteractors(neighborsOfBoth[i]);
      unsigned C_int_Ni_Count = nextClust.count();

      for (unsigned x = nextClust.find_first(); x < clusterGraph.numNodes(); x = nextClust.find_next(x)) {
	for (unsigned y = clusterGraph.getInteractors(x).find_first(); y < clusterGraph.numNodes();
	     y = clusterGraph.getInteractors(x).find_next(y)) {
	  if (( y != neighborsOfBoth[i] ) && (startClust[y] == 0)) {
	    ++T[y];
	    if ((y < neighborsOfBoth[i]) && (T[y] == C_int_Ni_Count)) {
	      return;
	    }
	  }
	}
      }

      // Calculate S[y] = |N(y) intersection with (C - N(i))| for y not in C
      vector<unsigned> S(clusterGraph.numNodes(), 0);
      boost::dynamic_bitset<unsigned long> C_not_Ni = startClust;
      C_not_Ni -= clusterGraph.getInteractors(neighborsOfBoth[i]);
      unsigned C_not_Ni_count = C_not_Ni.count();
      for (unsigned x = C_not_Ni.find_first(); x < clusterGraph.numNodes(); x = C_not_Ni.find_next(x)) {
	for (unsigned y = clusterGraph.getInteractors(x).find_first(); y < clusterGraph.numNodes();
	     y = clusterGraph.getInteractors(x).find_next(y)) {
	  if ((startClust[y] == 0) && (neighborsOfBothBits[y] == 1)) {
	    ++S[y];
	  }
	}
      }

      unsigned k = 0;
      unsigned last_jk = 0;
      for (unsigned jk = C_not_Ni.find_first(); jk < clusterGraph.numNodes(); jk = C_not_Ni.find_next(jk)) {
	boost::dynamic_bitset<unsigned long> Njk_not_C = clusterGraph.getInteractors(jk);
	for (unsigned y = Njk_not_C.find_first(); y < neighborsOfBoth[i]; y = Njk_not_C.find_next(y)) {
	  if ((T[y] == C_int_Ni_Count) && (neighborsOfBothBits[y] == 1)) {
	    if (y >= jk) {
	      --S[y];

	      // Might need to look here
	    } else if (((S[y] + k) == C_not_Ni_count) && (y >= last_jk)) {
	      return;
	    }
	  }
	}
	last_jk = jk;
	++k;
      }
      
      unsigned jp = last_jk;
      if (jp < (neighborsOfBoth[i]-1)) {
	return;
      }
      unsigned inNext = nextClust.find_first();
      for (unsigned y = clusterGraph.getInteractors(inNext).find_first(); y < clusterGraph.numNodes();
	   y = clusterGraph.getInteractors(inNext).find_next(y)) {
	if ((y < neighborsOfBoth[i]) && (startClust[y] == 0) && (T[y] == C_int_Ni_Count) && (S[y] == 0) && (jp < y)) {
	  return;
	}
      }
      nextClust[neighborsOfBoth[i]] = 1;
      updateCliques(clusterGraph, newClustersToAdd, nextClust, neighborsOfBoth, neighborsOfBothBits, i + 1, nodeIDsToNames);
    }
  }   

  // Return true if there is a cluster which contains this one.  Return false otherwise
  bool findClustsWithSeed(boost::dynamic_bitset<unsigned long> seedClust, graph_undirected_bitset & clusterGraph, 
			  vector<boost::dynamic_bitset<unsigned long> > & newClustersToAdd, vector<string> & nodeIDsToNames) {
    //boost::dynamic_bitset<unsigned long> newClust(clusterGraph.numNodes());
    //newClust[seedEdge.first] = 1;
    //newClust[seedEdge.second] = 1;
    
    // Find all newClusters Containing the seed cluster
    // unsigned newClustsContain = 0;
    //boost::dynamic_bitset<unsigned long> inNewClustWithThisEdge(clusterGraph.numNodes());
    for (unsigned i = 0; i < newClustersToAdd.size(); ++i) {
      if (seedClust.is_subset_of(newClustersToAdd[i])) {
	//++newClustsContain;
	return true;
	//inNewClustWithThisEdge |= newClustersToAdd[i];
      }
    }

    //if (newClustsContain == 0) {
    // Make a vector containing all nodes which are neighbors of all nodes contained in the seed cluster (called "Both" but really "All")
    vector<unsigned> neighborsOfBoth;
    vector<char> neighborsOfBothBits(clusterGraph.numNodes(), 0);
    neighborsOfBoth.reserve(clusterGraph.numNodes());
    for (unsigned i = 0; i < clusterGraph.numNodes(); ++i) {
      if ((seedClust[i] == 0) && (seedClust.is_subset_of(clusterGraph.getInteractors(i)))) {
	neighborsOfBoth.push_back(i);
	neighborsOfBothBits[i] = 1;
      }
    }

    //cout << "# Edge: " << seedEdge.first << "\t" << seedEdge.second << "\t" << neighborsOfBoth.size() << endl;
    //cout << "# Checking..." << nodeIDsToNames[seedEdge.first] << "\t" << nodeIDsToNames[seedEdge.second] << "\t" << neighborsOfBoth.size() << endl;
    if (neighborsOfBoth.size() > 0) {
      updateCliques(clusterGraph, newClustersToAdd, seedClust, neighborsOfBoth, neighborsOfBothBits, 0, nodeIDsToNames);
      return true;
    }
    //cout << "# Edge: " << nodeIDsToNames[seedEdge.first] << "\t" << nodeIDsToNames[seedEdge.second] << endl;
    return false;
  }


  void updateClustersWithEdges(vector<pair<pair<unsigned, unsigned>, double> > & edgesToAdd,
			       currentClusterClassBitset & currentClusters, graph_undirected_bitset & clusterGraph,
			       nodeDistanceObject & nodeDistances, unsigned & lastCurrent, vector<string> & nodeIDsToNames) {

    //cout << "# Adding " << edgesToAdd.size() << " edges" << endl;
    //vector<vector<char> > edgesCoveredByNew(vector<char>(clusterGraph.numNodes(), vector<char>(clusterGraph.numNodes(), 0)));
    unsigned long edgesToAddCounter = 0;
    if (edgesToAdd.size() != 0) {
      currentClusters.resetAllUnexplained();
    }
    while (edgesToAddCounter < edgesToAdd.size()) {
      
      unsigned long thisRoundCounter = edgesToAddCounter;
      
      vector<char> affectedNodes(clusterGraph.numNodes(), 0);
      //for (; (edgesToAddCounter == 0) || ((edgesToAddCounter != edgesToAdd.size()) && (edgesToAddCounter % 200000) != 0) ; ++edgesToAddCounter) {
      for (; edgesToAddCounter != edgesToAdd.size() ; ++edgesToAddCounter) {
	clusterGraph.addEdge(edgesToAdd[edgesToAddCounter].first.first, edgesToAdd[edgesToAddCounter].first.second);
	affectedNodes[edgesToAdd[edgesToAddCounter].first.first] = 1;
	affectedNodes[edgesToAdd[edgesToAddCounter].first.second] = 1;
      }
      
      unsigned numAff = 0;
      for (vector<char>::iterator thisIt = affectedNodes.begin(); thisIt != affectedNodes.end(); ++thisIt) {
	if (*thisIt) {
	  ++numAff;
	}
      }
      //cout << "# Nodes affected = " << numAff << endl;
      
      // For each affected node, find all possibly affected clusters
      vector<char> affectedClusters(currentClusters.maxClusterID(), 0);
      for (unsigned i = 0; i < affectedNodes.size(); ++i) {
	if (affectedNodes[i]) {
	  for (vector<unsigned long>::iterator cNodesIt = currentClusters.clustersWithNodeBegin(i);
	       cNodesIt != currentClusters.clustersWithNodeEnd(i); ++cNodesIt) {
	    affectedClusters[*cNodesIt] = 1;
	  }
	}
      }

      numAff = 0;
      for (vector<char>::iterator thisIt = affectedClusters.begin(); thisIt != affectedClusters.end(); ++thisIt) {
	if (*thisIt) {
	  ++numAff;
	}
      }
      //cout << "# Clusters affected = " << numAff << endl;
      
      
      vector<boost::dynamic_bitset<unsigned long> > newClustersToAdd;
      newClustersToAdd.reserve(5000);

      unsigned deleted = 0;
      for (unsigned i = 0; i < affectedClusters.size(); ++i) {
	if (affectedClusters[i]) {
	  
	  if (findClustsWithSeed(currentClusters.getElements(i), clusterGraph, newClustersToAdd, nodeIDsToNames)) {
	    // Cluster was not maximal - delete it
	    currentClusters.deleteCluster(i, nodeIDsToNames, false);
	    ++deleted;
	  }
	}
      }
      //cout << "# Num deleted here: " << deleted << endl;

      for ( ; thisRoundCounter < edgesToAddCounter; ++thisRoundCounter) {
	boost::dynamic_bitset<unsigned long> seedClust(clusterGraph.numNodes());
	seedClust[edgesToAdd[thisRoundCounter].first.first] = 1;
	seedClust[edgesToAdd[thisRoundCounter].first.second] = 1;
	if (!findClustsWithSeed(seedClust, clusterGraph, newClustersToAdd, nodeIDsToNames)) {
	  // seedClust is a maximal clique of just two nodes
	  Utils::insertInOrder(newClustersToAdd, seedClust);
	}
      }
      //cout << "# Found " << newClustersToAdd.size() << " new clusters to add" << endl;

      /*unsigned deleted = 0;
      for (unsigned i = 0; i < affectedClusters.size(); ++i) {
	if (affectedClusters[i]) {
	 
	  if (findClustsWithSeed(currentClusters.getElements(i), clusterGraph, newClustersToAdd, nodeIDsToNames)) {
	    // Cluster was not maximal - delete it
	    currentClusters.deleteCluster(i, nodeIDsToNames, false);
	    ++deleted;
	  }
	}
      }
      cout << "# Num deleted here: " << deleted << endl;*/

      for (vector<boost::dynamic_bitset<unsigned long> >::iterator clustersToAddIt = newClustersToAdd.begin();
           clustersToAddIt != newClustersToAdd.end(); ++clustersToAddIt) {
        currentClusters.addCluster(*clustersToAddIt,nodeDistances, nodeIDsToNames);
      }

      //cout << "Clusters after edges added: " << currentClusters.numCurrentClusters() << endl;
      /*cout << "# Added all edges. ";
      time (&end);
      dif = difftime(end,start);
      cout << "Time elapsed: " << dif << " seconds" << endl;
      time(&start);*/
      
      if ((currentClusters.numCurrentClusters() > 4*lastCurrent) || ((currentClusters.numCurrentClusters() > lastCurrent) && ((currentClusters.numCurrentClusters() - lastCurrent) > 1000))) {
	//cout << "CHECKING VALIDITY" << endl;
	//cout << "Num edges: " << clusterGraph.numEdges() << " Before delete: " << currentClusters.numCurrentClusters();
	performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	lastCurrent = currentClusters.numCurrentClusters();
	//cout << " After delete: " << lastCurrent << endl;
	if (lastCurrent < 25) {
	  lastCurrent = 25;
	}
      }

      //cout << "Clusters after edges added: " << currentClusters.numCurrentClusters() << endl;
      /*cout << "# Validity check after combining clusters done ";
      time (&end);
      dif = difftime(end,start);
      cout << "Time elapsed: " << dif << " seconds" << endl;*/
      ++edgesToAddCounter;
    }
  }

  void combineClusters(vector<pair<pair<unsigned long, unsigned long>, double> > & clustersToCombine,
		       currentClusterClassBitset & currentClusters, graph_undirected_bitset & clusterGraph,
		       unsigned & lastCurrent, vector<string> & nodeIDsToNames,
		       nodeDistanceObject & nodeDistances) {

    time_t start, end;
    time (&start);

    vector<vector<char> > edgeWillBeAdded(nodeDistances.numNodes(), vector<char>(nodeDistances.numNodes(), false));
    vector<pair<pair<unsigned, unsigned>, double> > edgesToAdd;
    edgesToAdd.reserve(clusterGraph.numEdges());
    
    // Sort clustersToCombine in order of descending weight
    sort(clustersToCombine.begin(), clustersToCombine.end(), compClustersToCombine);

    //cout << "# Sorted. ";
    time (&end);
    double dif = difftime(end,start);
    //cout << "Time elapsed: " << dif << " seconds" << endl;
    time (&start);

    //cout << currentClusters.maxClusterID() << endl;
    //graph_undirected clustersToCombineGraph(currentClusters.maxClusterID());
    //vector<char> clustersToBeDeleted(currentClusters.maxClusterID(), 0);
    
    double weight = clustersToCombine.front().second;
    unsigned beginThisWeight = 0;
    for (vector<pair<pair<unsigned long, unsigned long>, double> >::iterator clustersToCombineIt = clustersToCombine.begin();
	 clustersToCombineIt != clustersToCombine.end(); ++clustersToCombineIt) {
    
      //cout << "ADDING CLUSTER EDGE " << clustersToCombineIt->first.first << "\t" << clustersToCombineIt->first.second << endl;
      //clustersToCombineGraph.addEdge(clustersToCombineIt->first.first, clustersToCombineIt->first.second);
      //clustersToBeDeleted[clustersToCombineIt->first.first] = 1;
      //clustersToBeDeleted[clustersToCombineIt->first.second] = 1;
  
      if (clustersToCombineIt->second != weight) {
	if (edgesToAdd.size() != beginThisWeight) {
	  sort(edgesToAdd.begin() + beginThisWeight, edgesToAdd.end(), compEdgesToAdd);
	  beginThisWeight = edgesToAdd.size();
	}
	weight = clustersToCombineIt->second;
      }
      //cout << "Combining " << cluster1 << " " << cluster2 << endl;

      /*cout << "Combining\t";
      printCluster(currentClusters.getElements(clustersToCombineIt->first.first), nodeIDsToNames);
      cout << endl;
      printCluster(currentClusters.getElements(clustersToCombineIt->first.second), nodeIDsToNames);
      cout << endl << "-----" << endl;*/


      boost::dynamic_bitset<unsigned long> onlyIn1 = currentClusters.getElements(clustersToCombineIt->first.first);
      boost::dynamic_bitset<unsigned long> onlyIn2 = currentClusters.getElements(clustersToCombineIt->first.second);
      onlyIn1 -= currentClusters.getElements(clustersToCombineIt->first.second);
      onlyIn2 -= currentClusters.getElements(clustersToCombineIt->first.first);
      
      //currentClusters.setCurWeight(clustersToCombineIt->second);

      // ADD EDGES BETWEEN GENES THAT ARE ONLY IN ONE OF THE CLUSTERS TO THE CLUSTER GRAPH
      for (unsigned i = onlyIn1.find_first(); i < onlyIn1.size(); i = onlyIn1.find_next(i)) {
	boost::dynamic_bitset<unsigned long> newEdgesFrom1To2 = onlyIn2;
	newEdgesFrom1To2 -= clusterGraph.getInteractors(i);
	for (unsigned j = newEdgesFrom1To2.find_first(); j < newEdgesFrom1To2.size(); j = newEdgesFrom1To2.find_next(j)) {
	  //unsigned firstNode = i;
	  //unsigned secondNode = j;
	  if (!edgeWillBeAdded[i][j]) {
	    edgesToAdd.push_back(make_pair(make_pair(i,j), weight));
	    edgeWillBeAdded[i][j] = true;
	    edgeWillBeAdded[j][i] = true;
	    nodeDistances.setDistance(i,j,weight);
	  }
	}
      }
    }
    
    //cout << "# Got all " << edgesToAdd.size() << " edges. ";
    time (&end);
    dif = difftime(end,start);
    //cout << "Time elapsed: " << dif << " seconds" << endl;
    time (&start);

    updateClustersWithEdges(edgesToAdd, currentClusters, clusterGraph, nodeDistances, lastCurrent, nodeIDsToNames);

    //cout << "# Added all edges. ";
    time (&end);
    dif = difftime(end,start);
    //cout << "Time elapsed: " << dif << " seconds" << endl;
    return;
  }

  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  // dagConstruct::addMissingEdges
  // Takes clusters and combines highly overlapping clusters
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////

  bool addMissingEdges(currentClusterClassBitset & currentClusters, graph_undirected_bitset & clusterGraph, double density, double threshold, 
		       unsigned & lastCurrent, vector<string> & nodeIDsToNames, graph_undirected_bitset & realEdges, 
		       nodeDistanceObject & nodeDistances) {
    // clustersToCombine is a pair.  The first element of the pair is itself a pair, which is the ids of the two clusters to combine.
    // The second element of the pair is the weight of the combined edges
    vector<pair<pair<unsigned long, unsigned long>, double> > clustersToCombine;

    unsigned long maxClusterID = currentClusters.maxClusterID();
    vector<bool> clustersChecked(maxClusterID, false);
    vector<unsigned long> clustersToRecheck;
    for (unsigned long i = 0; i < maxClusterID; ++i) {
      if ((currentClusters.numElements(i) != 0) && currentClusters.isNew(i) && 
	  //((currentClusters.getClusterWeight(i) - currentClusters.getCurWeight()) > threshold)) {
	  (currentClusters.getThresh(i) > currentClusters.getCurWeight())) {
	clustersChecked[i] = true;
	for (unsigned long j = 0; j < maxClusterID; ++j) {
	  if ((j != i) && (currentClusters.numElements(j) != 0) &&
	      !clustersChecked[j]) {
	    if (isMinNodeDegreeMet(i,j,currentClusters,clusterGraph,density)) {
	      double newEdgeWeight = currentClusters.getClusterWeight(i);
	      if (currentClusters.getClusterWeight(j) < newEdgeWeight) {
		newEdgeWeight = currentClusters.getClusterWeight(j);
	      }
	      clustersToCombine.push_back(make_pair(make_pair(i,j), newEdgeWeight));
	      Utils::insertInOrder(clustersToRecheck, j);
	    }
	  }
	}
      }
    }

    while (clustersToRecheck.size()) {
      vector<unsigned long> clustersToRecheckAgain;
      for (vector<unsigned long>::iterator it = clustersToRecheck.begin(); it != clustersToRecheck.end(); ++it) {
	clustersChecked[*it] = true;
	for (unsigned long j = 0; j < maxClusterID; ++j) {
	  if ((*it != j) && (currentClusters.numElements(j) != 0) &&
	      !clustersChecked[j]) {
	    if (isMinNodeDegreeMet(*it,j,currentClusters,clusterGraph,density)) {
	      double newEdgeWeight = currentClusters.getClusterWeight(*it);
	      if (currentClusters.getClusterWeight(j) < newEdgeWeight) {
		newEdgeWeight = currentClusters.getClusterWeight(j);
	      }
	      //clustersToCombine.push_back(make_pair(make_pair(currentClusters.getElements(*it),currentClusters.getElements(j)), newEdgeWeight));
	      clustersToCombine.push_back(make_pair(make_pair(*it,j), newEdgeWeight));
	      Utils::insertInOrder(clustersToRecheckAgain, j);
	    }
	  }
	}
      }
      clustersToRecheck = clustersToRecheckAgain;
    }

    //cout << "# Combining " << clustersToCombine.size() << " pairs of clusters" << endl;
    /*for (unsigned count = 0; count < clustersToCombine.size(); ++count) {
      cout << clustersToCombine[count].second << endl;
      }*/

    if (clustersToCombine.size() > 0) {
      double curWeight = currentClusters.getCurWeight();
      //currentClusters.combining();
      combineClusters(clustersToCombine, currentClusters, clusterGraph, lastCurrent, nodeIDsToNames, nodeDistances);
      //currentClusters.doneCombining();
      currentClusters.setCurWeight(curWeight);
      return true;
    }
    return false;
  }
  
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  /// Uses a bitset representation of cluster graph                    ///
  /// and clusters to get maximal cliques                              ///
  /// Uses a threshold to deal with imperfect case                     ///
  /// Returns number of top level nodes                                ///
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  unsigned getFuzzyClustersBitsetThreshold(nodeDistanceObject & nodeDistances, map<string,unsigned> & nodeNamesToIDs,
					   vector<validClusterBitset> & validClusters, double threshold = 0, double density = 1) {
    time_t start, end;
    time (&start);
    unsigned largestCluster = 0;
    unsigned lastCurrent = 10;
    unsigned numRealEdgesAdded = 0;

    vector<string> nodeIDsToNames(nodeNamesToIDs.size(), string(""));
    for (map<string,unsigned>::iterator it = nodeNamesToIDs.begin(); it != nodeNamesToIDs.end(); ++it) {
      nodeIDsToNames[it->second] = it->first;
    }
    unsigned numNodes = nodeDistances.numNodes();
    graph_undirected_bitset clusterGraph(numNodes);
    graph_undirected_bitset realEdges(numNodes);
    double dt = nodeDistances.sortedDistancesBegin()->second; // Current threshold distance
    currentClusterClassBitset currentClusters(numNodes, clusterGraph.get_num_blocks(),dt, threshold);
    //currentClusters.setCurWeight(dt);
    bool addAtEnd = true;
    sortedDistanceStruct::iterator distanceIt = nodeDistances.sortedDistancesBegin();
    //for (sortedDistanceStruct::iterator distanceIt = nodeDistances.sortedDistancesBegin();
    //distanceIt != nodeDistances.sortedDistancesEnd(); ++distanceIt) {
    unsigned totalEdges = numNodes*(numNodes-1) / 2;
    while ((clusterGraph.numEdges() != totalEdges) && (distanceIt != nodeDistances.sortedDistancesEnd()) && (distanceIt->second >= Corrector::max()*threshold)) {

      vector<pair<pair<unsigned, unsigned>, double> > edgesToAdd;
      edgesToAdd.reserve(2000000);
      double addUntil = currentClusters.getNextThresh();
      if ((dt - Corrector::max()*threshold) > addUntil) {
	addUntil = dt - Corrector::max()*threshold;
      }
      //cout << "# " << distanceIt->second << "\t" << addUntil << "\t" << Corrector::max()*threshold << endl;
      while ((distanceIt != nodeDistances.sortedDistancesEnd()) && (distanceIt->second >= addUntil) && (distanceIt->second >= Corrector::max()*threshold)) {
	unsigned firstNode = distanceIt->first.first;
	unsigned secondNode = distanceIt->first.second;
	realEdges.addEdge(firstNode, secondNode);
	++numRealEdgesAdded;
	if (!clusterGraph.isEdge(firstNode, secondNode)) {
	  edgesToAdd.push_back(make_pair(make_pair(firstNode, secondNode), distanceIt->second));
	}
	++distanceIt;
      }
      updateClustersWithEdges(edgesToAdd, currentClusters, clusterGraph, nodeDistances, lastCurrent, nodeIDsToNames);

      //cout << "# numNew: " << currentClusters.numNew() << endl;
      
      //if (distanceIt->second != dt) {
      double last_dt = dt;
      if (distanceIt != nodeDistances.sortedDistancesEnd()) {
	dt = distanceIt->second;
      } else {
	dt = 0;
      }
      currentClusters.setCurWeight(dt);

      // THINK ABOUT WHETHER THIS IS STILL CORRECT
      double maxThresh = currentClusters.getMaxThresh();
      //cout << "# Max thresh: " << maxThresh << "\t";
      //cout << currentClusters.numNew() << endl;
      if ((dt < Corrector::max()*threshold) && (maxThresh == 0)) {
	addAtEnd = false;
	break;
      }
      
      //cout << "# New level " << currentClusters.getNextThresh() << "\t" << dt << endl;
      //if ((currentClusters.getMaxNewWeight() - distanceIt->second) > threshold) {
      
      // If there is a cluster which is now passing the threshold, check to see if it is valid or needs to be combined
      //cout << "# CHECK FOR CHECK " << currentClusters.getMaxThresh() << "\t" << dt << endl;
      if (maxThresh >= dt) {
	unsigned long numClustersBeforeDelete = currentClusters.numCurrentClusters();
	
	vector<unsigned long> newClustersSorted;
	if (density < 1) {
	  performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	  
	  // IN HERE WE NEED TO CHECK FOR MISSING EDGES BY COMBINING CLUSTERS INTO DENSE CLUSTERS
	  //cout << "# Adding missing edges...checking " << currentClusters.numCurrentClusters() << " cliques" << endl;
	  //currentClusters.combining();
	  bool newEdgesAdded = addMissingEdges(currentClusters, clusterGraph, density, threshold, lastCurrent, nodeIDsToNames, realEdges, nodeDistances);
	  while (newEdgesAdded == true) {
	    //cout << "# Rechecking" << endl;
	    performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	    newEdgesAdded = addMissingEdges(currentClusters, clusterGraph, density, threshold, lastCurrent, nodeIDsToNames, realEdges, nodeDistances);
	  }
	  
	  // Add current clusters to valid clusters if we have reached a new threshold distance
	  //currentClusters.removeNonMaximalClusters(clusterGraph, nodeDistances);
	  currentClusters.sortNewClusters(newClustersSorted);
	  //currentClusters.doneCombining();
	} else {
	  //currentClusters.removeNonMaximalClusters(clusterGraph, nodeDistances,nodeIDsToNames);
	  currentClusters.prepareForValidityCheck(newClustersSorted);
	}
	
	vector<char> idsChecked(currentClusters.maxClusterID(), 0);
	for (vector<unsigned long>::iterator newClusterIt = newClustersSorted.begin();
	     newClusterIt != newClustersSorted.end(); ++newClusterIt) {
	  bool isNecessary = false;
	  bool checkForFinal = false;
	  double clustWeight = currentClusters.getClusterWeight(*newClusterIt);
	  //if (currentClusters.isNew(*newClusterIt) && ((clustWeight - distanceIt->second) > threshold)) {
	  if (currentClusters.isNew(*newClusterIt) && (currentClusters.getThresh(*newClusterIt) > distanceIt->second)) {
	    checkForFinal = true;
	    currentClusters.setCheckedFinal(*newClusterIt);
	  }
	  
	  if (currentClusters.checkClusterFinalValidity(*newClusterIt,isNecessary, idsChecked, checkForFinal)) {
	    //validClusters.push_back(validClusterBitset(currentClusters.getElements(*newClusterIt), 0, (1 - clustWeight) / 2.0));
	    /*cout << "# Adding valid cluster:\t" << *newClusterIt << "\t" << currentClusters.numElements(*newClusterIt) << "\t" << clustWeight << endl;
	      cout << "# History:\t";
	      for (vector<pair<double,unsigned> >::iterator weightIt = currentClusters.historicalWeightsBegin(*newClusterIt);
	      weightIt != currentClusters.historicalWeightsEnd(*newClusterIt); ++weightIt) {
	      cout << weightIt->second << "," << weightIt->first << "\t";
	      }
	      cout << endl;
	      currentClusters.clearHistoricalWeights(*newClusterIt);*/
	    
	    validClusters.push_back(validClusterBitset(currentClusters.getElements(*newClusterIt), 0, clustWeight));
	    //cout << "# Valid cluster:\t";
	    printCluster(currentClusters.getElements(*newClusterIt), nodeIDsToNames);
	    //currentClusters.setClusterValid(currentClusters.getElements(*newClusterIt));
	    currentClusters.setClusterValid(*newClusterIt);
	    //cout << "\t" << clustWeight << "\t" << currentClusters.getNumUniquelyUnexplainedEdges(*newClusterIt) << "\t" << currentClusters.getThresh(*newClusterIt) << "\t" << dt << endl;
	    if (validClusters.back().numElements() > largestCluster) {
	      largestCluster = validClusters.back().numElements();
	    }
	    currentClusters.setOld(*newClusterIt);
	    
	  } else if (!isNecessary) {
	    //bool alreadyCombining = currentClusters.isCombiningNow();
	    //currentClusters.combining();
	    //cout << "# Deleting here" << endl;
	    currentClusters.deleteCluster(*newClusterIt,nodeIDsToNames,false);
	    
	    /*if (!alreadyCombining) {
	      currentClusters.doneCombining();
	      }*/
	  } else if (isNecessary) { // && currentClusters.isNew(*newClusterIt) && !currentClusters.wasNecessary(*newClusterIt) && checkForFinal) {
	    currentClusters.setNecessary(*newClusterIt);
	    // MAY HAVE TO THINK ABOUT THIS MORE
	    //currentClusters.setNumExtensions(*newClusterIt, 0);
	    /*cout << "Found necessary:\t";
	      printCluster(currentClusters.getElements(*newClusterIt), nodeIDsToNames);
	      cout << "\t" << clustWeight << endl;*/
	  }
	}
	lastCurrent = currentClusters.numCurrentClusters();
	if (lastCurrent < 25) {
	  lastCurrent = 25;
	}
	
	/*
	cout << "# dt: " << last_dt << endl;
	cout << "# Next dt: " << dt << endl;
	cout << "# Num current clusters before delete: " << numClustersBeforeDelete << endl;
	cout << "# Num current clusters: " << currentClusters.numCurrentClusters() << endl;
	cout << "# Num valid clusters: " << validClusters.size() << endl;
	cout << "# Largest cluster: " << largestCluster << endl;
	cout << "# Num edges in clusterGraph: " << clusterGraph.numEdges() << endl;
	cout << "# Num real edges: " << numRealEdgesAdded << endl;
	cout << "# Num edges inferred: " << clusterGraph.numEdges() - numRealEdgesAdded << endl;
	* */
	time (&end);
	double dif = difftime(end,start);
	//cout << "# Time elapsed: " << dif << " seconds" << endl;
      }
    }

    /*if ((currentClusters.numCurrentClusters() > 4*lastCurrent) || ((currentClusters.numCurrentClusters() > lastCurrent) && ((currentClusters.numCurrentClusters() - lastCurrent) > 1000))) {
	//cout << "CHECKING VALIDITY" << endl;
	performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	lastCurrent = currentClusters.numCurrentClusters();
	if (lastCurrent < 25) {
	  lastCurrent = 25;
	}
	//cout << "LAST CURRENT: " << lastCurrent << endl;
	}*/

    /*vector<pair<pair<unsigned, unsigned>, double> > edgesToAdd;
      edgesToAdd.reserve(2000000);
      double addUntil = currentClusters.getNextThresh();
      if ((dt - Corrector::max()*threshold) > addUntil) {
	addUntil = dt - Corrector::max()*threshold;
      }
      //cout << "# " << distanceIt->second << "\t" << addUntil << endl;
      while (distanceIt->second >= addUntil) {
	unsigned firstNode = distanceIt->first.first;
	unsigned secondNode = distanceIt->first.second;
	realEdges.addEdge(firstNode, secondNode);
	++numRealEdgesAdded;
	if (!clusterGraph.isEdge(firstNode, secondNode)) {
	  //addEdgeAndUpdateClusters(firstNode, secondNode, clusterGraph, currentClusters, nodeIDsToNames, nodeDistances);
	  edgesToAdd.push_back(make_pair(make_pair(firstNode, secondNode), distanceIt->second));
	}
	++distanceIt;
      }
      --distanceIt;
      updateClustersWithEdges(edgesToAdd, currentClusters, clusterGraph, nodeDistances, lastCurrent, nodeIDsToNames);
      */
    
    // Add clusters at end to valid clusters
    if (addAtEnd) {
      //cout << "# Adding at end" << endl;
      double last_dt = dt;
      dt = 0;
      currentClusters.setCurWeight(0);
      vector<unsigned long> newClustersSorted;
      if (density < 1) {
	performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	
	// IN HERE WE NEED TO CHECK FOR MISSING EDGES BY COMBINING CLUSTERS INTO DENSE CLUSTERS
	//cout << "# Adding missing edges...checking " << currentClusters.numCurrentClusters() << " clusters" << endl;
	//currentClusters.combining();
	bool newEdgesAdded = addMissingEdges(currentClusters, clusterGraph, density, threshold, lastCurrent, nodeIDsToNames, realEdges, nodeDistances);
	while (newEdgesAdded == true) {
	  //cout << "# Rechecking" << endl;
	  performValidityCheck(currentClusters, clusterGraph, nodeDistances, nodeIDsToNames);
	  newEdgesAdded = addMissingEdges(currentClusters, clusterGraph, density, threshold, lastCurrent, nodeIDsToNames, realEdges, nodeDistances);
	}
	//currentClusters.removeNonMaximalClusters(clusterGraph, nodeDistances);
	currentClusters.sortNewClusters(newClustersSorted);
	//currentClusters.doneCombining();
      } else {
	//currentClusters.removeNonMaximalClusters(clusterGraph, nodeDistances, nodeIDsToNames);
	currentClusters.prepareForValidityCheck(newClustersSorted);
      }
      
      vector<char> idsChecked(currentClusters.maxClusterID(), 0);
      for (vector<unsigned long>::iterator newClusterIt = newClustersSorted.begin();
	   newClusterIt != newClustersSorted.end(); ++newClusterIt) {
	bool isNecessary = false;
	bool checkForFinal = false;
	double clustWeight = currentClusters.getClusterWeight(*newClusterIt);
	//if (currentClusters.isNew(*newClusterIt) && (clustWeight > threshold)) {
	if (currentClusters.isNew(*newClusterIt) && (currentClusters.getThresh(*newClusterIt) > 0)) {
	  checkForFinal = true;
	}
	
	if (checkForFinal && currentClusters.checkClusterFinalValidity(*newClusterIt,isNecessary, idsChecked, checkForFinal)) {
	  
	  /*cout << "# Adding valid cluster:\t" << *newClusterIt << "\t" << currentClusters.numElements(*newClusterIt) << "\t" << clustWeight << endl;
	  cout << "# History:\t";
	  for (vector<pair<double,unsigned> >::iterator weightIt = currentClusters.historicalWeightsBegin(*newClusterIt);
	       weightIt != currentClusters.historicalWeightsEnd(*newClusterIt); ++weightIt) {
	    cout << weightIt->second << "," << weightIt->first << "\t";
	  }
	  cout << endl;
	  currentClusters.clearHistoricalWeights(*newClusterIt);*/

	  validClusters.push_back(validClusterBitset(currentClusters.getElements(*newClusterIt), 0, clustWeight));
	  //cout << "# Valid cluster:\t";
	  printCluster(currentClusters.getElements(*newClusterIt), nodeIDsToNames);
	  //currentClusters.setClusterValid(currentClusters.getElements(*newClusterIt));
	  currentClusters.setClusterValid(*newClusterIt);
	  //cout << "\t" << clustWeight << "\t" << currentClusters.getNumUniquelyUnexplainedEdges(*newClusterIt) << "\t" << currentClusters.getThresh(*newClusterIt) << "\t" << dt << endl;
	  if (validClusters.back().numElements() > largestCluster) {
	    largestCluster = validClusters.back().numElements();
	  }
	  currentClusters.setOld(*newClusterIt);

	} else if (checkForFinal && !isNecessary) {
	  currentClusters.deleteCluster(*newClusterIt, nodeIDsToNames, false);
	}
      }

	/*
      cout << "# dt: " << last_dt << endl;
      cout << "# Num current clusters: " << currentClusters.numCurrentClusters() << endl;
      cout << "# Num valid clusters: " << validClusters.size() << endl;
      cout << "# Largest cluster: " << largestCluster << endl;
      cout << "# Num edges in clusterGraph: " << clusterGraph.numEdges() << endl;
      */
    }
    return currentClusters.numCurrentClusters();
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
  // Take list of valid clusters and turn it into an ontology, assuming the perfect case
  // where all edges in a cluster are identical
  //////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////

  void clustersToDAG(vector<validClusterBitset> validClusters, DAGraph & ontology, unsigned numTerminalNodes) {
    //unsigned numTerminalNodes = nodeDistances.numNodes();
    
    sort(validClusters.begin(), validClusters.end());
    

    ontology.reserveNodes(ontology.numNodes() + validClusters.size() + 1);
    unsigned clustCount = 0;
    vector<validClusterBitset>::iterator clustersIt = validClusters.begin();
    unsigned firstNodeID = ontology.addNode();
    clustersIt->setID(firstNodeID);
    for (unsigned i = 0; i < numTerminalNodes; ++i) {
      if (clustersIt->isElement(i)) {
	ontology.addEdge(firstNodeID, i);
      }
    }
    ontology.setWeight(clustersIt->getWeight(),firstNodeID);
    double geneWeight = clustersIt->getWeight();
    ++clustersIt;
    ++clustCount;
    for ( ; clustersIt != validClusters.end(); ++clustersIt) {
      //cout << clustCount << endl;
      unsigned newNodeID = ontology.addNode();
      clustersIt->setID(newNodeID);

      boost::dynamic_bitset<unsigned long> unaccountedFor = clustersIt->getElements();

      for (vector<validClusterBitset>::reverse_iterator possibleDescendentsIt(clustersIt); 
	   possibleDescendentsIt != validClusters.rend(); ++possibleDescendentsIt) {
	if (possibleDescendentsIt->numElements() < clustersIt->numElements()) {
	  // A cluster is a child if it is not already a 
	  // descendent of the current cluster via some other cluster,
	  // and all of its genes are contained in the possible ancestor cluster,
	  
	  bool isAlreadyDescendent = ontology.isDescendent(possibleDescendentsIt->getID(), newNodeID);
	  if ((!isAlreadyDescendent) &&
	      (isClusterAncestor(possibleDescendentsIt->getElements(), clustersIt->getElements(), unaccountedFor))) {
	    //cout << "adding edge " << newNodeID << "\t" << possibleDescendentsIt->getID() << endl;
	    ontology.addEdge(newNodeID, possibleDescendentsIt->getID());
	  }
	}
      }
      
      // Add edges to genes which are not contained in any of the child nodes
      for (unsigned i = 0; i < numTerminalNodes; ++i) {
	if (unaccountedFor[i] == 1) {
	  ontology.addEdge(newNodeID, i);
	}
      }
      
      // Set weight of cluster
      ontology.setWeight(clustersIt->getWeight(),newNodeID);
      if (geneWeight < clustersIt->getWeight()) {
	geneWeight = clustersIt->getWeight();
      }
      //cout << ontology.getName(newNodeID) << "\t" << ontology.getWeight(newNodeID) << endl;
      ++clustCount;
    }

    if (geneWeight > 1) {
      ontology.setGeneWeight(geneWeight);
    } else {
      ontology.setGeneWeight(1);
    }

    //cout << "Finished creating DAG - adding root node" << endl;
    unsigned firstTopLevelNode;
    bool firstTopLevelFound = false;
    bool secondTopLevelFound = false;
    long rootID = -1;
    
    unsigned numTopLevel = 0;
    for (vector<DAGNode>::iterator nodeIt = ontology.nodesBegin(); nodeIt != ontology.nodesEnd(); ++nodeIt) {
      if ((nodeIt->numParents() == 0) && (nodeIt->getID() != rootID)) {
	//cout << nodeIt->getID() << " " << nodeIt->getName() << " is top level" << endl;
	++numTopLevel;
	//cout << numTopLevel << endl;
	if (firstTopLevelFound == false) {
	  firstTopLevelFound = true;
	  firstTopLevelNode = nodeIt->getID();
	} else if (secondTopLevelFound == false) {
	  secondTopLevelFound = true;
	  unsigned curItPos = nodeIt->getID();

	  rootID = ontology.addNode();
	  // Adding node may invalidate nodeIt, so fix it
	  nodeIt = ontology.nodesBegin();
	  nodeIt += curItPos;

	  ontology.setWeight(0, rootID);
	  ontology.addEdge(rootID, firstTopLevelNode);
	  ontology.addEdge(rootID, nodeIt->getID());
	} else {
	  if (rootID != nodeIt->getID()) {
	    ontology.addEdge(rootID, nodeIt->getID());
	  }
	}
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  // dagConstruct::constructDAG(...)
  // Main clustering function - cluster input graph into Ontology
  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  void constructDAG(graph_undirected & input_graph, DAGraph & ontology, nodeDistanceObject & nodeDistances, double threshold, double density) {
    //ontology = DAGraph();
    map<string,unsigned> geneNamesToIDs;

    // Add all nodes in the input network to the ontology as genes
    for (vector<Node>::iterator nodesIt = input_graph.nodesBegin(); nodesIt != input_graph.nodesEnd(); ++nodesIt) {
      geneNamesToIDs[nodesIt->getName()] = nodesIt->getID();
      ontology.addNode(nodesIt->getName(), geneNamesToIDs);
    }

    nodeDistances = nodeDistanceObject(input_graph);

    vector<validClusterBitset> validClusters;
    dagConstruct::getFuzzyClustersBitsetThreshold(nodeDistances, geneNamesToIDs, validClusters, threshold, density);

    // Got the clusters - build ontology assuming any clusters wholly contained in others are descendents
    //cout << "# Num clusters: " << validClusters.size() << endl;
    dagConstruct::clustersToDAG(validClusters, ontology, geneNamesToIDs.size());
    return;
  }
}
#endif //DAG_CONSTRUCT
