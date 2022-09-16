# tests for ACTIONet network construction
# make sure you import the data if you don't have it locally: ./get_test_datasets.sh

import scanpy as sc

from _ACTIONet import reduce_kernel
from ACTIONet.main import run_ACTIONet

# read the data
adata = sc.read_h5ad("CellSIUS_Wegmann19/normalized_input.h5ad")

# perform reduction
reduce_kernel(adata)

# build network
print("JSD/KNN")
an_jsd_knn = run_ACTIONet(adata, network_metric="jsd", network_algorithm="knn", copy=True)
an_jsd_knn.write_h5ad("clecius.jsd.knn.h5ad")
print("JSD/K*NN")
an_jsd_kstarnn = run_ACTIONet(adata, network_metric="jsd", network_algorithm="k*nn", copy=True)
an_jsd_kstarnn.write_h5ad("celsius.jsd.kstarnn.h5ad")
print("L2/K*NN")
an_l2_kstarnn = run_ACTIONet(adata, network_metric="l2", network_algorithm="k*nn", copy=True)
an_l2_kstarnn.write_h5ad("celsius.l2.kstarnn.h5ad")
print("IP/K*NN")
an_ip_kstarnn = run_ACTIONet(adata, network_metric="ip", network_algorithm="k*nn", copy=True)
an_ip_kstarnn.write_h5ad("celsius.ip.kstarnn.h5ad")
