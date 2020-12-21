from typing import Optional

from anndata import AnnData

from . import network as nt
from . import preprocessing as pp


def run_ACTIONet(
    adata: AnnData,
    k_max: Optional[int] = 30,
    min_cells_per_archetype: Optional[int] = 2,
    min_specificity_z_threshold: Optional[int] = -3,
    network_density: Optional[int] = 1,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 500,
    layout_in_parallel: Optional[bool] = True,
    n_threads: Optional[int] = 0,
    unification_alpha: Optional[float] =  0.99, 
    unification_outlier_threshold: Optional[float] = 2, 
    unification_sim_threshold: Optional[float] = 0, 
    footprint_alpha: Optional[float] = 0.85,
    max_iter_ACTION: Optional[int] = 50,
    full_trace: Optional[bool] = False,
    copy: Optional[bool] = False,
):
    adata = adata.copy() if copy else adata

    # Run ACTION
    C, H = pp.ACTION(
        adata,
        k_min=2,
        k_max=k_max,
        n_threads=n_threads,
        max_it=max_iter_ACTION,
        min_delta=1e-300,
    )

    # Prune nonspecific and/or unreliable archetypes
    pp.prune_archetypes(
        adata,
        C,
        H,
        min_specificity_z_threshold=min_specificity_z_threshold,
        min_cells=min_cells_per_archetype,
    )

    # Build ACTIONet
    nt.build_network(
        adata,
        density=network_density,
        n_threads=n_threads,
        mutual_edges_only=mutual_edges_only,
    )

    # Layout ACTIONet
    if not layout_in_parallel:
        nt.layout_network(
            adata,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=1,
        )
    else:
        nt.layout_network(
            adata,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=n_threads,
        )

    # Identiy equivalent classes of archetypes and group them together
    pp.unify_archetypes(
        adata,
        unification_alpha, 
        unification_outlier_threshold, 
        unification_sim_threshold, 
        n_threads)

    # Use graph core of global and induced subgraphs to infer centrality/quality of
    # each cell
    nt.compute_archetype_core_centrality(adata)

    # Smooth archetype footprints
    nt.compute_network_diffusion(adata, alpha=footprint_alpha, n_threads=n_threads)

    # Compute gene specificity for each archetype
    pp.compute_archetype_feature_specificity(adata)

    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=round(layout_epochs / 5),
        n_threads=1,
    )

    return adata if copy else None
