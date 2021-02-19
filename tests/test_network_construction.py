#tests for ACTIONet network construction
#make sure you import the data if you don't have it locally: ./get_test_datasets.sh

import ACTIONet as an
import scanpy as sc
 
#read the data
adata=sc.read_h5ad("CellSIUS_Wegmann19/normalized_input.h5ad")

#perform reduction
an.pp.reduce_kernel(adata)

#build network
an.run_ACTIONet(adata)



