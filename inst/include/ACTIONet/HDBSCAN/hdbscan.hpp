#pragma once
#include <ACTIONet.h>
#include <string>
#include <vector>
#include "hdbscanParameters.hpp"
#include "hdbscanResult.hpp"
#include "hdbscanRunner.hpp"
#include "outlierScore.hpp"

using namespace std;

class Hdbscan

{
 private:
  hdbscanResult result;

 public:
  Hdbscan(mat& X);

  void import_arma(mat& X);

  vector<vector<double> > dataset;

  std::vector<int> labels_;

  std::vector<int> normalizedLabels_;

  std::vector<outlierScore> outlierScores_;

  std::vector<double> membershipProbabilities_;

  uint32_t noisyPoints_;

  uint32_t numClusters_;

  void execute(int minPoints, int minClusterSize, string distanceMetric);

  void displayResult();
};
