from typing import Optional

import numpy as np
from anndata import AnnData

from . import annotation as annot
from . import decomposition as decomp
from . import network as net


def run_ACTIONet(
        adata: AnnData,
        k_min: Optional[int] = 2,
        k_max: Optional[int] = 30,
        layer_key: Optional[str] = None,
        reduction_key: Optional[str] = "ACTION",
        net_key_out: Optional[str] = "ACTIONet",
        min_cells_per_archetype: Optional[int] = 2,
        max_iter_ACTION: Optional[int] = 50,
        specificity_th: Optional[int] = -3,
        network_metric: Optional[str] = "jsd",
        network_algorithm: Optional[str] = "k*nn",
        network_density: Optional[int] = 1,
        network_k: Optional[int] = 30,
        mutual_edges_only: Optional[bool] = True,
        layout_compactness: Optional[int] = 50,
        layout_epochs: Optional[int] = 1000,
        layout_algorithm: Optional[str] = "tumap",
        layout_in_parallel: Optional[bool] = True,
        unification_th: Optional[float] = 0,
        footprint_key: Optional[str] = "archetype_footprint",
        footprint_alpha: Optional[float] = 0.85,
        thread_no: Optional[int] = 0,
        seed: Optional[int] = 0,
        copy: Optional[bool] = False,
        ):
    adata = adata.copy() if copy else adata

    if layer_key is not None:
        if layer_key not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_key + "']. ")
        S = adata.layers[layer_key]
    else:
        S = adata.X

    S = S.astype(dtype=np.float64)

    if reduction_key not in adata.obsm.keys():
        raise ValueError("Did not find adata.obsm['" + reduction_key + "']. ")

    alg_name = layout_algorithm.lower()
    if alg_name not in ["umap", "tumap"]:
        raise ValueError("'layout_algorithm' must be 'tumap' or 'umap'.")

    decomp.runACTIONMR(
            data=adata,
            k_min=k_min,
            k_max=k_max,
            max_iter=max_iter_ACTION,
            min_delta=1e-300,
            specificity_th=specificity_th,
            min_cells_per_archetype=min_cells_per_archetype,
            unification_th=unification_th,
            thread_no=thread_no,
            reduction_key=reduction_key,
            return_W=False,
            return_raw=False,
            copy=False,
            )

    # Build ACTIONet
    net.build(
            data=adata,
            algorithm=network_algorithm,
            distance_metric=network_metric,
            density=network_density,
            mutual_edges_only=mutual_edges_only,
            k=network_k,
            net_key_out=net_key_out,
            H_key="H_stacked",
            thread_no=thread_no,
            copy=False,
            return_raw=False,
            )

    # Uses `reduction_key` for initialization of coordinates
    net.layout(
            data=adata,
            algorithm=alg_name,
            initial_coordinates=None,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            thread_no=thread_no if layout_in_parallel else 1,
            reduction_key=reduction_key,
            net_key=net_key_out,
            seed=seed,
            copy=False,
            return_raw=False,
            )

    # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
    net.centrality(
            adata,
            algorithm="localized_coreness",
            annotations_key="assigned_archetype",
            copy=False,
            return_raw=False,
            )

    # Smooth archetype footprints
    net.diffusion(
            data=adata,
            scores_key="H_stacked",
            smoothed_scores_key=footprint_key,
            alpha_val=footprint_alpha,
            thread_no=thread_no,
            copy=False,
            return_raw=False,
            )

    # Compute gene specificity for each archetype
    annot.compute_archetype_feature_specificity(
            adata=adata,
            S=S,
            H=None,
            layer_key=None,
            footprint_key=footprint_key,
            copy=False,
            return_raw=False,
            )

    # net.construct_backbone(
    #     adata=adata,
    #     footprint=None,
    #     net_key=net_key_out,
    #     footprint_key=footprint_key,
    #     scale=True,
    #     layout_compactness=layout_compactness,
    #     layout_epochs=round(layout_epochs / 5),
    #     alpha_val=footprint_alpha,
    #     thread_no=1,
    #     seed=seed,
    #     copy=False
    # )

    return adata if copy else None
