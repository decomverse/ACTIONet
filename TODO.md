# TODO

## Primary
* Update and test geneset scoring
* Update and test marker selection (one vs all, pair-wise)
* Test, debug, and run code on Windows - submit to Bioconductor
* Visualization: Interactive cell selection
* Visualization: "gating" on archetypes
* Test batch-corrected ACTION, batch-orthogonalization, and other methods fo BC
* Package markers as a JSON file(s)
* Package pathways/genesets as JSON file(s)
* Simplify and streamline installation
* Bookdown tutorial(s)
* Add back backbone + spanner + visualization
* Test unification
* Speed-up gene/cell filtering
* Fix `correct.cell.annotations` to take character

## To be tested
* `transform_layout()` for embedding datasets into ACTIONet plot
* `compute_pseudo_bulk_per_archetype, compute_pseudo_bulk_per_cluster, compute_pseudo_bulk_per_ind` functions for pseudo-bulk.
* `compute_AA_coreset()` for AA coreset construction + wAA for fast sketching
* `sgd2_layout_weighted()` for [S_GD2 layout](https://github.com/jxz12/s_gd2) layout
* `compute_sparse_network_diffusion` for  **L1-regularized PageRank** algorithm algorithm from: ["Variational perspective on local graph clustering"](https://github.com/kfoynt/LocalGraphClustering)
* Variance-adjusted Limma + Add weighted mean/std (first version -- Sebastian)
* weighted AA
* Test and debug Online AA


## Extentions
* subACTIONet using weighted AA
* reACTION
* Supervized AA
* Construct cell-state ontology graph using ([circle packing](http://jeromefroe.github.io/circlepackeR/))
* https://github.com/astamm/fdahotelling


## Done
* Deflated SVD with A*B being prior archetypal analysis (done, used in batch ortho)
* Backbone (archetype) network construction algorithm using PageRank+ACTIONet (in unification)
* Write, test, and release docker-lite and docker-full - also install Rstudio Servver/Jupyter (done). 


