# Installation
## Setting up environment
First decide on which compiler/BLAS/LAPACK configuration to use. Depending on the choice, create (if it doesn't exist) a .R/Makevars (linux and mac) or .R/Makevars.win (windows) file in your home directory and copy the content of the following files into it:

* **Linux/Windows (OS) + oneAPI (dpcpp compiler + MKL)**: `src/Makevars_lib/Makevars.dpcpp` **--** Best / Fastest (if using Intel architechture)
* **Linux/Mac/Windows (OS) + gcc (compiler) + MKL**: `src/Makevars_lib/Makevars.gcc_MKL` **--** Second best / Fast enough
* **Mac (OS) + clang (compiler) + accelerate**: `src/Makevars_lib/Makevars.clang` **--** Easy on Mac / Fast-ish
* **Linux/Mac/Windows (OS) + gcc (compiler) + default BLAS/LAPACK**: `src/Makevars_lib/Makevars.gcc` **--** Easiest / [likely] Slowest

#### > PS:
* For further installation instruction for oneAPI and/or oneMKL, visit: [Installation Guide for Intel® oneAPI Toolkits](https://software.intel.com/content/www/us/en/develop/articles/installation-guide-for-intel-oneapi-toolkits.html).
* The fastest speed will be achieved if R itself is also compiled with oneMKL/oneAPI. For compiling R, follow instructions in [Install R with MKL](X)
* If your Mac architecture doesn't use Intel processors, using is not suggested.


## Installing ACTIONet R package
If insalling using oneMKL/oneAPI, make sure correct paths are set by running `/opt/intel/inteloneapi/setvars.sh`.

### From within R

```r
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet")
```

### Directly from the clone

```bash
git clone git@github.com:shmohammadi86/ACTIONet.git
cd ACTIONet
R CMD INSTALL .
```


## Install optional packages
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
ACTIONet interfaces to scran normalization via `scran.normalize()` function.

* [Linnorm](https://bioconductor.riken.jp/packages/3.4/bioc/html/Linnorm.html) is another commonly used normalization technique. You can install it via:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

biocLite("Linnorm")

```


ACTIONet interfaces to scran normalization via `linnorm.normalize()` function.


#### General Single-cell processing
* Another common package used in R for single-cell analysis is [Seurat](https://satijalab.org/seurat/). You can install the stable version from CRAN:

```r
install.packages('Seurat')

```
or the latest developmental version from github:

```r
install.packages('devtools')
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
```
ACTIONet contains a method to import data from a Seurat object: `import.ace.from.Seurat()`.

# Running ACTIONet
Here is a simple example to get you started:

```r
# Download example dataset from the 10X Genomics website
download.file('http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5', 'pbmc_10k_v3.h5') 

require(ACTIONet)
# Run ACTIONet
ace = import.ace.from.10X.h5('pbmc_10k_v3.h5', prefilter = T, min_cells_per_feat = 50, min_umis_per_cell = 1000)
ace = reduce.ace(ace)
ACTIONet_results = run.ACTIONet(ace)
ace = ACTIONet_results$ace

# Annotate cell-types
data("curatedMarkers_human")
markers = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
annot.out = annotate.cells.from.archetypes.using.markers(ace, markers)
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





# Additional tutorials
You can access ACTIONet tutorials from:
1. [ACTIONet framework at a glance (human PBMC 3k dataset)](http://compbio.mit.edu/ACTIONet/min_intro.html)
2. [Introduction to the ACTIONet framework (human PBMC Granja et al. dataset)](http://compbio.mit.edu/ACTIONet/intro.html)
3. [Introduction to cluster-centric analysis using the ACTIONet framework](http://compbio.mit.edu/ACTIONet/clustering.html)
4. [To batch correct or not to batch correct, that is the question!](http://compbio.mit.edu/ACTIONet/batch.html)
5. [PortingData: Import/export options in the ACTIONet framework](http://compbio.mit.edu/ACTIONet/porting_data.html)
6. [Interactive visualization, annotation, and exploration](http://compbio.mit.edu/ACTIONet/annotation.html)
7. [Constructing cell-type/cell-state-specific networks using SCINET](http://compbio.mit.edu/ACTIONet/scinet.html)

You can also find a [step-by-step guide](http://compbio.mit.edu/ACTIONet/guide.html) to learning the core functionalities of the ACTIONet framework, as well as a detailed description of the [ACTIONetExperiment (ACE)](http://compbio.mit.edu/ACTIONet/ace.html) object.
