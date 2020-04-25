#pragma once
#include <string>
#include <vector>
#include "hdbscanRunner.hpp"
#include "hdbscanParameters.hpp"
#include "hdbscanResult.hpp"
#include "outlierScore.hpp"
#include <ACTIONet.h>

using namespace std;


class Hdbscan

{

private:

	hdbscanResult result;

public:
	Hdbscan(mat& X);

	void import_arma(mat& X);	
	
	vector < vector <double > > dataset;

	std::vector<int> labels_;

	std::vector<int> normalizedLabels_;

	std::vector<outlierScore>outlierScores_;

	std::vector <double> membershipProbabilities_;

	uint32_t noisyPoints_;

	uint32_t numClusters_;


	void execute(int minPoints, int minClusterSize, string distanceMetric);

	void displayResult();

	

};

