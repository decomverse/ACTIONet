#ifndef GRAPH_UNDIRECTED_BITSET
#define GRAPH_UNDIRECTED_BITSET

#include "boost.h"

class graph_undirected_bitset {
 public:
  graph_undirected_bitset(unsigned numNodes) {
    graph = vector<boost::dynamic_bitset<unsigned long> >(numNodes, boost::dynamic_bitset<unsigned long>(numNodes));
    numTotalEdges = 0;
  }

  inline unsigned get_bits_per_block() {
    return graph[0].bits_per_block;
  }

  inline unsigned get_num_blocks() {
    return graph[0].num_blocks();
  }

  inline void addEdge(unsigned node1, unsigned node2) {
    graph[node1][node2] = 1;
    graph[node2][node1] = 1;
    ++numTotalEdges;
  }

  inline unsigned numEdges() {
    return numTotalEdges;
  }

  const boost::dynamic_bitset<> & getInteractors(unsigned nodeID) {
    return graph[nodeID];
  }

  inline bool isEdge(unsigned node1, unsigned node2) {
    return graph[node1].test(node2);
  }

  static void getIntersection(boost::dynamic_bitset<unsigned long> & set1, boost::dynamic_bitset<unsigned long> & set2, boost::dynamic_bitset<unsigned long> & result) {
    result = set1;
    result &= set2;
  }

  inline unsigned numNodes() {
    return graph.size();
  }

  void printAll(vector<string> & nodeIDsToNames) {
    cout << "@ PRINTING GRAPH" << endl;
    unsigned numEdgesFound = 0;
    for (unsigned i = 0; i < graph.size(); ++i) {
      unsigned long j = graph[i].find_first();
      while (j < i) {
	j = graph[i].find_next(j);
      }
      if (j != graph[i].npos) {
	cout << "@ " << nodeIDsToNames[i] << "\t" << nodeIDsToNames[j];
	while ((j = graph[i].find_next(j)) < graph[i].npos) {
	  cout << "," << nodeIDsToNames[j];
	  ++numEdgesFound;
	}
	cout << endl;
      }
      if (numEdgesFound == numTotalEdges) {
	break;
      }
    }
    cout << "@ DONE PRINTING GRAPH" << endl;
  }

 private:
  vector<boost::dynamic_bitset<> > graph;
  unsigned numTotalEdges;
};


#endif // GRAPH_UNDIRECTED_BITSET

