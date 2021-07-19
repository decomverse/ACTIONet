#include "EuclideanDistance.hpp"
#include <cmath>
#include <cstdint>
#include <vector>
double EuclideanDistance::computeDistance(std::vector<double> attributesOne,
                                          std::vector<double> attributesTwo) {
  double distance = 0;
  for (uint32_t i = 0; i < attributesOne.size() && i < attributesTwo.size();
       i++) {
    distance += ((attributesOne[i] - attributesTwo[i]) *
                 (attributesOne[i] - attributesTwo[i]));
  }

  return sqrt(distance);
}
