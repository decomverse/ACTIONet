#Installation

### Setting Up the Environment (Preinstallation)
**For Linux Users**
For the optimal performance on Intel-based architectures, installing [Intel Math Kernel Library (MKL)](https://software.intel.com/content/www/us/en/develop/articles/intel-math-kernel-library-intel-mkl-2020-install-guide.html) is **highly** recommended. After installing, make sure `MKLROOT` is defined by running the [setvars](https://software.intel.com/content/www/us/en/develop/documentation/using-configuration-file-for-setvars-sh/top.html) script.

**Install library dependencies**
To install the `ACTIONet` dependencies on debian-based linux machines, run:

```bash
sudo apt install libhdf5-dev libsuitesparse-dev libcurl4-openssl-dev libssl-dev libxml2-dev
```

For Mac-based systems use [brew](https://brew.sh/):

```bash
brew install hdf5 suite-sparse
```
Some of the dependent R packages also require [XQuartz](https://www.xquartz.org) to be installed.

### Installing ACTIONet R Package

#### Using `devtools`
This is the easiest way to install the package, and it automatically installs all package dependencies (except optional packages described below):

```r
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-devel")

```

#### Directly from the clone
This is the more flexible approach and allows for easily pulling the latest updates/changes. This can be done as follows:

* **Instaling dependencies**: Unlike devtools, `R` script does not automatically install dependencies. To install all dependencies, run:

```r
install.packages(c('Matrix', 'Rcpp', 'RcppArmadillo', 'R.utils', 'hdf5r', 'plotly', 'ggpubr', 'corrplot', 'wordcloud', 'threejs', 'RColorBrewer'))
BiocManager::install(c("SingleCellExperiment", "ComplexHeatmap"))
```

* **Clone `ACTIONet` repository**:
If you don't already have git, install it first.

On (Debian) **Linux-based  machines**, run:

```bash
sudo apt-get install git
```

For **Mac-based machines**, run:

```bash
brew install git
```

Now clone a fresh copy of the repository:

```bash
git clone https://github.com/shmohammadi86/ACTIONet.git
```


* **Install `ACTIONet`**:
ACTIONet contains many different branchs. For installing the stable version, switch to the `R-release` branch:

```bash
cd ACTIONet
git checkout R-release
```

If you want to install the latest updates, switch to the `R-devel` branch:

```bash
cd ACTIONet
git checkout R-devel
```

and now you can install ACTIONet using the following command in the `ACTIONet` directory:

```bash
R CMD INSTALL .
```

### Install optional packages
#### Batch correction
* [batchelor](https://bioconductor.org/packages/release/bioc/html/batchelor.html): This Implements a variety of methods for batch correction of single-cell (RNA sequencing) data, including mutually-nearest neighbor (MNN) method. You can install it using bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("batchelor")
```

ACTIONet includes interface to run MNN: `reduce.and.batch.correct.sce.fastMNN()`.


* [Harmony](https://github.com/immunogenomics/harmony): Harmony is another popular batch-correction method that has direct interface implemented in the ACTIONet framework:

```r
install.packages("devtools")
devtools::install_github("immunogenomics/harmony")
```

ACTIONet includes interface to run harmony: `reduce.and.batch.correct.sce.Harmony()`.


#### Normalization & QC
* [scater](http://bioconductor.org/packages/release/bioc/html/scater.html)/[scran](https://bioconductor.org/packages/release/bioc/html/scran.html) packages provide great set of tools for normalization and quality-control of single-cell datasets stored as a `SingleCellExperiment` format. You can instal them using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("scater", "scran")
```
ACTIONet interfaces to scran normalization via `normalize.scran()` function.

* [Linnorm](https://bioconductor.riken.jp/packages/3.4/bioc/html/Linnorm.html) is another commonly used normalization technique. You can install it via:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

biocLite("Linnorm")

```


ACTIONet interfaces to scran normalization via `normalize.Linnorm()` function.



# Running ACTIONet
**Note** If you are using `MKL`, make sure to properly [set the number of threads](https://software.intel.com/content/www/us/en/develop/documentation/mkl-macos-developer-guide/top/managing-performance-and-memory/improving-performance-with-threading/techniques-to-set-the-number-of-threads.html) used prior to running `ACTIONet`.

## Example Run
Here is a simple example to get you started:

```r
# Download example dataset from the 10X Genomics website
download.file('http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5', 'pbmc_10k_v3.h5')

require(ACTIONet)
# Run ACTIONet
ace = import.ace.from.10X.h5('pbmc_10k_v3.h5', prefilter = T, min_cells_per_feat = 0.01, min_feats_per_cell = 1000)
ace = normalize.ace(ace)
ace = reduce.ace(ace)
ace = run.ACTIONet(ace)

# Annotate cell-types
data("curatedMarkers_human")
markers = curatedMarkers_human$Blood$PBMC$Ding2019$marker.genes
annot.out = annotate.cells.using.markers(ace, markers)
ace$celltypes = annot.out$Labels

# Visualize output
plot.ACTIONet(ace, "celltypes", transparency.attr = ace$node_centrality)

# Export results as AnnData
ACE2AnnData(ace, fname = "pbmc_10k_v3.h5ad")
```
## Visualizing results using cellxgene

ACTIONet framework introduces an extension of the `SingleCellExperiment` object that can be closely mapped to [AnnData](https://anndata.readthedocs.io/en/stable/index.html) object. In fact, output of ACTIONet in the python implementation is internally stored as as `AnnData` object, and R `ACE` objects can be imported from/exported to `AnnData` using functions `AnnData2ACE()` and `ACE2AnnData()` functions, respectively. `AnnData` objects can be directly loaded into [cellxgene](https://github.com/chanzuckerberg/cellxgene) package, an open-source viewer for interactive single-cell data visualization. `cellxgene` can be installed as:

```bash
pip install cellxgene

```

Then to visualize the results of ACTIONet, run:
```bash
cellxgene launch pbmc_10k_v3.h5ad
```

where *pbmc_10k_v3.h5ad* is the name of the file we exported using `ACE2AnnData()` function.


# Running `ACTIONet` with Docker

With Docker, no installation is required and users can directly start experimenting with `ACTIONet` instantly! Sounds easy, right? Here is how it goes.

* Download [Docker](https://store.docker.com/search?offering=community&type=edition) (if you don't already have, which in most linux distributions you do!)

* Run `ACTIONet` docker:

```bash
$ docker run -e USER=<USER> -e PASSWORD=<PASSWORD> -p 8787:8787 actionet/actionet:mro
```
Replace `<USER>` and `<PASSWORD>` with a username and password of your choice.

If you wish to access your local data inside the docker, modify the command as:

```bash
$ docker run -v /your/data/file/path/:/data -w /data -e USER=<USER> -e PASSWORD=<PASSWORD> -p 8787:8787 actionet/actionet:mini
```

* Connect to the `RStudio Server` through your desktops browser on `http://127.0.0.1:8787`, and then enter "rstudio" as the username and password when prompted. `ACTIONet` and all of its dependencies, as well as a few other packages for single-cell analysis,  are already installed in this environment for you.

* Have a cup of coffee and enjoy running ACTIONet!


# Additional tutorials
You can access ACTIONet tutorials from:
1. [ACTIONet framework at a glance (human PBMC 3k dataset)](http://compbio.mit.edu/ACTIONet/tutorials/mini_intro.html)
2. [Introduction to the ACTIONet framework (human PBMC Granja et al. dataset)](http://compbio.mit.edu/ACTIONet/tutorials/intro.html)
3. [Introduction to cluster-centric analysis using the ACTIONet framework](http://compbio.mit.edu/ACTIONet/tutorials/clustering.html)
4. [To batch correct or not to batch correct, that is the question!](http://compbio.mit.edu/ACTIONet/tutorials/batch.html)
5. [PortingData: Import/export options in the ACTIONet framework](http://compbio.mit.edu/ACTIONet/tutorials/porting_data.html)
6. [Interactive visualization, annotation, and exploration](http://compbio.mit.edu/ACTIONet/tutorials/annotation.html)
7. [Constructing cell-type/cell-state-specific networks using SCINET](http://compbio.mit.edu/ACTIONet/tutorials/scinet.html)

You can also find a [step-by-step guide](http://compbio.mit.edu/ACTIONet/tutorials/guide.html) to learning the core functionalities of the ACTIONet framework.
