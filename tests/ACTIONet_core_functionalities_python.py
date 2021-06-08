from ACTIONet import *

# Imports PBMC 3k dataset from scanpy
import scanpy as sc
pbmc = sc.datasets.pbmc3k_processed()
ACE = pbmc
ACE

# Performs kernel reduction
reduce_kernel(ACE)

# Runs ACTION for a range of archetypes (2..30 by default)
C, H = run_ACTION(ACE)

# Preprocesses archetypes to prune likely noisy archetypes.
prune_archetypes(ACE, C, H)

# Uses k*-NN together with JS metric on archetypes to construct cell-cell network (ACTIONet graph)
build_ACTIONet(ACE, mutual_edges_only = True)

# Computes 2D/3D embedding of cells, as well as their de novo coloring
layout_ACTIONet(ACE)

# Identifies and unifies redundant archetypes into equivalent classes and assigns cells to each class
unify_archetypes(ACE)

# Compute centrality of each vertex within the subgraph induced by assigned cells to each archetype
compute_archetype_core_centrality(ACE)


sc.plotting.embedding(ACE, "ACTIONet2D", color="louvain")
sc.plotting.embedding(ACE, "ACTIONet2D", color="assigned_archetypes")


# Evaluate clustering quality w.r.t. Leiden clusters (not the best, but just a sanity check)
celltypes = pbmc.obs.louvain
from sklearn.metrics.cluster import normalized_mutual_info_score
normalized_mutual_info_score(celltypes, ACE.obs["assigned_archetypes"])

ACE.write("/home/shahin/Dropbox/Projects/SingleCell/repositories/ACTION/meta_repo_python/tests/ACTIONet_out_python.h5ad")
