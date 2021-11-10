#!/bin/bash
apt-get update \
  && apt-get install -y --no-install-recommends libhdf5-dev libsuitesparse-dev libnss3 xvfb

Rscript -e "devtools::install_github('shmohammadi86/ACTIONetExperiment', ask = FALSE, upgrade = 'always')"
Rscript -e "devtools::install_github('shmohammadi86/ACTIONet', ref = 'R-devel', ask = FALSE, upgrade = 'always')"
