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

# Additional tutorials
You can access ACTIONet tutorials from:
* X


# Visualizing ACTIONet results using cellxgene
ACTIONet framework introduces an extension of the `SingleCellExperiment` object that can be closely mapped to [AnnData](https://anndata.readthedocs.io/en/stable/index.html) object. In fact, output of ACTIONet in the python implementation is internally stored as as `AnnData` object, and R `ACE` objects can be imported from/exported to `AnnData` using functions `AnnData2ACE()` and `ACE2AnnData()` functions, respectively. `AnnData` objects can be directly loaded into [cellxgene](https://github.com/chanzuckerberg/cellxgene) package, an open-source viewer for interactive single-cell data visualization. `cellxgene` can be installed as:

```bash
pip install cellxgene

```

and `*.h5ad` files exported via `ACE2AnnData()` can be loaded into `cellxgene` as:

```r
cellxgene launch [ACTIONet.h5ad]

```

where `ACTIONet.h5ad` can be replaced with the filename of interest. Furthermore, `AnnData` objects are native to [Scanpy](https://scanpy.readthedocs.io/en/latest/index.html) package. To install `Scanpy` and its dependencies, run:

```bash
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leiden

pip install scanpy
```

and use the following to import and have an initial inspection of the results:

```python
import scanpy as sc

ACTIONet = sc.read("ACTIONet.h5ad")
sc.pl.embedding(ACTIONet, "ACTIONet2D", color="assigned_archetype")
```

The key point to remember is that ACTIONet framework is inherently designed to avoid disjoint clustering, but the above listed command visualizes the discretized form of archetypes for illustration purposes. 