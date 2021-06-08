import pdb
import ACTIONet as an
import scanpy as sc
#read the data
adata=sc.read_h5ad("fetal_liver_100k+/ACTIONet.h5ad")

#perform reduction
an.pp.reduce_kernel(adata)
