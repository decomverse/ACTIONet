import scanpy as sc

import ACTIONet as an

# read the data
adata = sc.read_h5ad("fetal_liver_100k+/ACTIONet.h5ad")

# perform reduction
an.pp.reduce_kernel(adata)
