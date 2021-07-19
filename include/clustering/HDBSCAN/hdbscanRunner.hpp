#pragma once
#include "hdbscanParameters.hpp"
#include "hdbscanResult.hpp"
class hdbscanRunner {
 public:
  static hdbscanResult run(hdbscanParameters parameters);
};
