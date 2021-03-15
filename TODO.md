# TODO

## Primary
* Test, debug, and run code on Windows
* Write, test, and release docker-lite and docker-full - also install Rstudio Servver/Jupyter. 1- Minimal (Alpine)
* Write, test, and release singularity-lite and singularity-full
* "gating" on archetypes
* Interactive cell selection


## New features to be tested
* `renormalize_input_matrix` for post-ACTIONet normalization
* `compute_sparse_network_diffusion` for  **L1-regularized PageRank** algorithm algorithm from: ["Variational perspective on local graph clustering"](https://github.com/kfoynt/LocalGraphClustering)
* `transform_layout()` for embedding datasets into ACTIONet plot
* `compute_pseudo_bulk_per_archetype, compute_pseudo_bulk_per_cluster, compute_pseudo_bulk_per_ind` functions for pseudo-bulk.
* `sgd2_layout_weighted()` for [S_GD2 layout](https://github.com/jxz12/s_gd2) layout
* `HDBSCAN.clustering` for clustering ACTIONet
* `compute_AA_coreset()` for AA coreset construction + wAA for fast sketching
* Harmonized action to be checked

## Extensions
### Tier 1
* Variance-adjusted Limma + Add weighted mean/std
* Deflated SVD with A*B being prior archetypal analysis
* subACTIONet using weighted AA
* Backbone (archetype) network construction algorithm using PageRank+ACTIONet
* reACTION
* Autocorrelation to find continuous cell-states + joint analysis with multi-resolution application of Leiden to identify clustered archetypes
* Construct cell-state ontology graph using ([circle packing](http://jeromefroe.github.io/circlepackeR/))
* Test and debug Online AA

## Tier 2
* Implement the **cluster refinement** algorithm from: ["A simple and strongly-local flow-based method for cut improvement"](https://github.com/kfoynt/LocalGraphClustering)
* Incorporate CellMESH (https://github.com/shunfumao/cellmesh)
* Multi-resolution batch-correction inspired by Harmony
* Supervized marker suggestion

## Minor
* Fix annotate with single cell type bug
* Merge filter functions
