FROM bioconductor/bioconductor_docker
ARG GITHUB_PAT

LABEL name="actionet/r_release" \
	version=1.0.0 \
	vendor="ACTIONet" \
	maintainer="shahin.mohammadi@gmail.com" \
	description="ACTIONet single-cell framework." \
	license="GPLv2"


ENV DEBIAN_FRONTEND=noninteractive \
	APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=1 \
	R_VERSION=$R_VERSION

## Install ACTIONet dependencies
RUN apt-get update \
	&& apt-get install -y --no-install-recommends libhdf5-dev libsuitesparse-dev libnss3 xvfb libglpk40

## Install package dependencies
RUN Rscript -e 'pkg_dep<-c("batchelor", "BiocParallel", "ComplexHeatmap", "ggplot2", "ggpubr", "ggrepel", "gplots", "grid", "hdf5r", "Matrix", "matrixStats", "pals", "plotly", "RColorBrewer", "Rcpp", "RcppArmadillo", "seriation", "SingleCellExperiment"); requireNamespace("BiocManager"); BiocManager::install(pkg_dep,ask=F);'

## Install ACTIONet
RUN Rscript -e "remotes::install_github('shmohammadi86/ACTIONetExperiment', ask = FALSE, upgrade = 'never'); remotes::install_github('shmohammadi86/ACTIONet', ref = 'R-release', ask = FALSE, upgrade = 'never');"

CMD ["/init"]
