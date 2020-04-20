import scanpy as sc
import numpy as np
from ACTIONet import *
from sklearn.metrics.cluster import normalized_mutual_info_score

pbmc = sc.datasets.pbmc3k_processed()
labels = pbmc.obs.louvain
S = np.transpose(pbmc.X)

reduction_out = reduce_kernel(S)
S_r = reduction_out['S_r']
S_r.shape

ACTION_out = run_ACTION(S_r, k_max = 10)

pruning_out = prune_archetypes(C_trace = ACTION_out["C"], H_trace = ACTION_out["H"])

G = build_ACTIONet(pruning_out["H_stacked"])
G.shape

vis_out = layout_ACTIONet(G, S_r)


unification_out = unify_archetypes(G, S_r, pruning_out["C_stacked"], pruning_out["H_stacked"])


specificity_out = compute_feature_specificity(S, unification_out["H_unified"])

