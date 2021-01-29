import numpy as np
from typing import Optional

from anndata import AnnData

from . import network as nt
from . import preprocessing as pp


def run_ACTIONet(
        adata: AnnData,
        k_max: Optional[int] = 30,
        layer_name: Optional[str] = None,
        reduction_name: Optional[str] = "ACTION",
        net_name_out: Optional[str] = "ACTIONet",
        min_cells_per_archetype: Optional[int] = 2,
        max_iter_ACTION: Optional[int] = 50,
        min_specificity_z_threshold: Optional[int] = -3,
        network_density: Optional[int] = 1,
        mutual_edges_only: Optional[bool] = True,
        layout_compactness: Optional[int] = 50,
        layout_epochs: Optional[int] = 1000,
        layout_algorithm: Optional[int] = 0,
        layout_in_parallel: Optional[bool] = True,
        unification_alpha: Optional[float] = 0.99,
        unification_outlier_threshold: Optional[float] = 2,
        unification_sim_threshold: Optional[float] = 0,
        footprint_alpha: Optional[float] = 0.85,
        thread_no: Optional[int] = 0,
        seed: Optional[int] = 0,
        copy: Optional[bool] = False,
):
    adata = adata.copy() if copy else adata

    if layer_name is not None:
        if layer_name not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_name + "']. ")
        S = adata.layers[layer_name]
    else:
        S = adata.X

    S = S.astype(dtype=np.float64)

    if reduction_name not in adata.obsm.keys():
        raise ValueError("Did not find adata.obsm['" + reduction_name + "']. ")
    else:
        S_r = adata.obsm[reduction_name].astype(dtype=np.float64)

    # Run ACTION
    ACTION_out = pp.ACTION(
        data=S_r,  # data must be `n_obs` × `n_vars/red_dim`
        reduction_name=None,
        k_min=2,
        k_max=k_max,
        thread_no=thread_no,
        max_it=max_iter_ACTION,
        min_delta=1e-300,
    )

    # Prune nonspecific and/or unreliable archetypes
    pp.prune_archetypes(
        ACTION_out["C"],
        ACTION_out["H"],
        adata=adata,
        min_specificity_z_threshold=min_specificity_z_threshold,
        min_cells=min_cells_per_archetype,
        copy=False
    )

    # Build ACTIONet
    nt.build_network(
        adata,
        net_name_out=net_name_out,
        density=network_density,
        thread_no=thread_no,
        mutual_edges_only=mutual_edges_only,
        copy=False,
        return_raw=False
    )
    # adata.obsp[net_name_out] = G

    nt.layout_network(
        adata=adata,
        G=None,
        S_r=None,
        reduction_name=reduction_name,
        net_name=net_name_out,
        compactness_level=layout_compactness,
        n_epochs=layout_epochs,
        layout_alg=layout_algorithm,
        thread_no=thread_no if layout_in_parallel else 1,
        seed=seed,
        copy=False,
        return_raw=False
    )

    # Identiy equivalent classes of archetypes and group them together
    pp.unify_archetypes(
        adata=adata,
        G=None,
        S_r=None,
        C_stacked=None,
        alpha_val=unification_alpha,
        outlier_threshold=unification_outlier_threshold,
        sim_threshold=unification_sim_threshold,
        thread_no=thread_no,
        copy=False,
        return_raw=False
    )

    # Use graph core of global and induced subgraphs to infer centrality/quality of
    # each cell
    nt.compute_archetype_core_centrality(adata)

    # Smooth archetype footprints
    nt.compute_network_diffusion(adata, alpha=footprint_alpha, thread_no=thread_no)

    # Compute gene specificity for each archetype
    pp.compute_archetype_feature_specificity(adata)

    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=round(layout_epochs / 5),
        thread_no=1,
    )

    return adata if copy else None
