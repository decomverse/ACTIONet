#include "hdbscan.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

Hdbscan::Hdbscan(mat& X) { import_arma(X); }

// Each row is a data point
void Hdbscan::import_arma(mat& X) {
  vector<vector<double> > dataset;

  vector<double> row(X.n_cols);
  for (int i = 0; i < X.n_rows; i++) {
    for (int j = 0; j < X.n_cols; j++) {
      row[j] = X(i, j);
    }
    dataset.push_back(row);
  }

  this->dataset = dataset;

  return;
}

void Hdbscan::execute(int minPoints, int minClusterSize,
                      string distanceMetric) {
  // Call The Runner Class here
  hdbscanRunner runner;
  hdbscanParameters parameters;
  uint32_t noisyPoints = 0;
  set<int> numClustersSet;
  map<int, int> clustersMap;
  vector<int> normalizedLabels;

  parameters.dataset = this->dataset;
  parameters.minPoints = minPoints;
  parameters.minClusterSize = minClusterSize;
  parameters.distanceFunction = distanceMetric;
  this->result = runner.run(parameters);
  this->labels_ = result.labels;
  this->outlierScores_ = result.outliersScores;
  for (uint32_t i = 0; i < result.labels.size(); i++) {
    if (result.labels[i] == 0) {
      noisyPoints++;
    } else {
      numClustersSet.insert(result.labels[i]);
    }
  }
  this->numClusters_ = numClustersSet.size();
  this->noisyPoints_ = noisyPoints;
  int iNdex = 1;
  for (auto it = numClustersSet.begin(); it != numClustersSet.end(); it++) {
    clustersMap[*it] = iNdex++;
  }
  for (int i = 0; i < labels_.size(); i++) {
    if (labels_[i] != 0)
      normalizedLabels.push_back(clustersMap[labels_[i]]);
    else if (labels_[i] == 0) {
      normalizedLabels.push_back(-1);
    }
  }
  this->normalizedLabels_ = normalizedLabels;
  this->membershipProbabilities_ = result.membershipProbabilities;
}

// void Hdbscan::displayResult() {
//   hdbscanResult result = this->result;
//   uint32_t numClusters = 0;
//
//   cout << "HDBSCAN clustering for " << this->dataset.size() << " objects."
//        << endl;
//
//   for (uint32_t i = 0; i < result.labels.size(); i++) {
//     cout << result.labels[i] << " ";
//   }
//
//   cout << endl << endl;
//
//   cout << "The Clustering contains " << this->numClusters_ << " clusters with "
//        << this->noisyPoints_ << " noise Points." << endl;
// }
