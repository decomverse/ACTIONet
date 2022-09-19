import scanpy as sc

from _ACTIONet import reduce_kernel

# read the data
adata = sc.read_h5ad("fetal_liver_100k+/ACTIONet.h5ad")

# perform reduction
reduce_kernel(adata)
