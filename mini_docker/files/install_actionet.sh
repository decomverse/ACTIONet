#!/bin/bash
Rscript -e "devtools::install_github('shmohammadi86/ACTIONetExperiment', ask = FALSE, upgrade = 'never')"
Rscript -e "devtools::install_github('shmohammadi86/ACTIONet', ref = 'R-devel', ask = FALSE, upgrade = 'never')"
