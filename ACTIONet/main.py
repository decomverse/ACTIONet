from typing import Optional

import numpy as np
from anndata import AnnData

import _ACTIONet as _an
from ACTIONet import postprocessing as po
from ACTIONet.decomposition.actionmr import runACTIONMR
from ACTIONet.network.build import build
from ACTIONet.network.centrality import centrality
from ACTIONet.network.diffusion import diffusion
from ACTIONet.network.layout import layout
from ACTIONet.tools.utils_public import normalize_reduction, scale_matrix


def run_ACTIONet(
    adata: AnnData,
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    layer_key: Optional[str] = None,
    reduction_key: Optional[str] = None,
    normalization: Optional[int] = 1,
    net_key_out: Optional[str] = "ACTIONet",
    min_cells_per_archetype: Optional[int] = 2,
    max_iter_ACTION: Optional[int] = 50,
    specificity_th: Optional[int] = -3,
    network_metric: Optional[str] = "jsd",
    network_algorithm: Optional[str] = "k*nn",
    network_density: Optional[int] = 1,
    network_k: Optional[int] = 30,
    mutual_edges_only: Optional[bool] = True,
    layout_epochs: Optional[int] = 100,
    layout_algorithm: Optional[str] = "umap",
    layout_presmooth_network: Optional[bool] = False,
    layout_sim2dist: Optional[int] = 2,
    layout_min_dist: Optional[float] = 1.0,
    layout_spread: Optional[float] = 1.0,
    layout_learning_rate: Optional[float] = 1.0,
    layout_in_parallel: Optional[bool] = True,
    unification_backbone_density: Optional[float] = 0.5,
    unification_resolution: Optional[float] = 1.0,
    unification_min_cluster_size: Optional[int] = 3,
    footprint_key: Optional[str] = "archetype_footprint",
    footprint_alpha: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    seed: Optional[int] = 0,
    backbbone_density: Optional[float] = 0.5,
    copy: Optional[bool] = False,
):
    adata = adata.copy() if copy else adata

    if layer_key is None and "default_assay" in adata.uns["metadata"].keys():
        layer_key = adata.uns["metadata"]["default_assay"]

    if layer_key is not None:
        if layer_key not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_key + "']. ")
        S = adata.layers[layer_key]
    else:
        S = adata.X
    S = S.astype(dtype=np.float64)

    if reduction_key is None:
        if "default_reduction" in adata.uns["metadata"].keys():
            reduction_key = adata.uns["metadata"]["default_reduction"]
        else:
            print("reduction_slot is None and 'default_reduction' does  not exist in adata.uns['metadata']. Reseting to `ACTION` as default")
            reduction_key = "ACTION"

    if reduction_key not in adata.obsm.keys():
        raise ValueError("Did not find adata.obsm['" + reduction_key + "']. ")

    normalize_reduction(adata=adata, reduction_key=reduction_key, normalization=normalization)

    reduction_key = reduction_key + "_normalized"

    alg_name = str(layout_algorithm).lower()
    if alg_name not in ["umap", "tumap"]:
        raise ValueError("'layout_algorithm' must be 'tumap' or 'umap'.")

    runACTIONMR(
        data=adata,
        k_min=int(str(k_min)),
        k_max=int(str(k_max)),
        max_iter=max_iter_ACTION,
        min_delta=1e-300,
        specificity_th=specificity_th,
        min_cells_per_archetype=min_cells_per_archetype,
        unification_backbone_density=unification_backbone_density,
        unification_min_cluster_size=unification_min_cluster_size,
        unification_resolution=unification_resolution,
        thread_no=thread_no,
        reduction_key=reduction_key,
        return_W=False,
        return_raw=False,
        copy=False,
    )

    # Build ACTIONet
    build(
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

    # Smooth archetype footprints
    diffusion(
        data=adata,
        scores_key="H_unified",
        smoothed_scores_key=footprint_key,
        alpha_val=footprint_alpha,
        thread_no=thread_no,
        copy=False,
        return_raw=False,
    )

    # Uses `reduction_key` for initialization of coordinates
    svd_out = _an.IRLB_SVD_full(
        scale_matrix(adata.obsm["archetype_footprint"]),
        dim=3,
        iters=1000,
        seed=0,
        verbose=0,
    )
    initial_coordinates = scale_matrix(svd_out["u"][:, 0:3])

    adata.obsm["H_stacked"]
    layout(
        data=adata,
        initial_coordinates=initial_coordinates,
        algorithm=alg_name,
        presmooth_network=layout_presmooth_network,
        sim2dist=layout_sim2dist,
        spread=layout_spread,
        min_dist=layout_min_dist,
        learning_rate=layout_learning_rate,
        n_epochs=layout_epochs,
        thread_no=thread_no if layout_in_parallel else 1,
        reduction_key=reduction_key,
        net_key=net_key_out,
        seed=seed,
        copy=False,
        return_raw=False,
    )

    # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
    centrality(
        adata,
        algorithm="localized_coreness",
        labels="assigned_archetype",
        copy=False,
        return_raw=False,
    )

    # Compute gene specificity for each archetype
    po.archetypes.feature_specificity(
        adata=adata,
        S=S,
        H=None,
        layer_key=None,
        footprint_key=footprint_key,
        copy=False,
        return_raw=False,
    )

    # Construct archetype graph (backbone)
    po.archetypes.construct_backbone(adata, density=backbbone_density)

    return adata if copy else None


def rerun_layout(
    adata: AnnData,
    layout_epochs: Optional[int] = 100,
    layout_algorithm: Optional[str] = "umap",
    layout_presmooth_network: Optional[bool] = False,
    layout_sim2dist: Optional[int] = 2,
    layout_min_dist: Optional[float] = 1.0,
    layout_spread: Optional[float] = 1.0,
    layout_learning_rate: Optional[float] = 1.0,
    net_key: Optional[str] = "ACTIONet",
    footprint_key: Optional[str] = "archetype_footprint",
    thread_no: Optional[int] = 0,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False,
):
    adata = adata.copy() if copy else adata

    # Uses `reduction_key` for initialization of coordinates
    svd_out = np.linalg.svd(scale_matrix(adata.obsm[footprint_key]).T)
    initial_coordinates = scale_matrix(svd_out[2][:, 0:3])

    adata.obsm["H_stacked"]
    layout(
        data=adata,
        initial_coordinates=initial_coordinates,
        algorithm=layout_algorithm,
        presmooth_network=layout_presmooth_network,
        sim2dist=layout_sim2dist,
        spread=layout_spread,
        min_dist=layout_min_dist,
        learning_rate=layout_learning_rate,
        n_epochs=layout_epochs,
        thread_no=thread_no,
        net_key=net_key,
        seed=seed,
        copy=False,
        return_raw=False,
    )

    return adata if copy else None
