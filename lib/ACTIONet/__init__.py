__author__ = ', '.join([
    'Shahin Mohammadi'
])
__email__ = ', '.join([
    'shahin.mohammadi@gmail.com'
])


# Load API
from ._decomposition import reduce_kernel, run_simplex_regression, run_SPA, run_AA, run_ACTION, prune_archetypes, unify_archetypes
from ._network_tools import build_ACTIONet, layout_ACTIONet, compute_archetype_core_centrality
from ._statistics import compute_archetype_feature_specificity, compute_cluster_feature_specificity


from anndata import AnnData
