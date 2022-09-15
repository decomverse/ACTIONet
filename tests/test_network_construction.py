# tests for ACTIONet network construction
# make sure you import the data if you don't have it locally: ./get_test_datasets.sh

import scanpy as sc

import ACTIONet as an

# read the data
adata = sc.read_h5ad("CellSIUS_Wegmann19/normalized_input.h5ad")

# perform reduction
an.pp.reduce_kernel(adata)

# build network
print("JSD/KNN")
an_jsd_knn = an.run_ACTIONet(adata, distance_metric="jsd", nn_approach="knn", copy=True, k=10)
an_jsd_knn.write_h5ad("clecius.jsd.knn.h5ad")
print("JSD/K*NN")
an_jsd_kstarnn = an.run_ACTIONet(adata, distance_metric="jsd", nn_approach="k*nn", copy=True)
an_jsd_kstarnn.write_h5ad("celsius.jsd.kstarnn.h5ad")
print("L2/K*NN")
an_l2_kstarnn = an.run_ACTIONet(adata, distance_metric="l2", nn_approach="k*nn", copy=True)
an_ls_kstarnn.write_h5ad("celsius.l2.kstarnn.h5ad")
print("IP/K*NN")
an_ip_kstarnn = an.run_ACTIONet(adata, distance_metric="ip", nn_approach="k*nn", copy=True)
an_ip_kstarnn.write_h5ad("celsius.ip.kstarnn.h5ad")
