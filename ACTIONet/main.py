import numpy as np
import pandas as pd
from typing import Optional
from natsort import natsorted

from anndata import AnnData

from . import network as net
from . import preprocessing as pp
from . import decomposition as decomp


def run_ACTIONet(
    adata: AnnData,
    k_max: Optional[int] = 30,
    layer_key: Optional[str] = None,
    reduction_key: Optional[str] = "ACTION",
    net_key: Optional[str] = "ACTIONet",
    min_cells_per_archetype: Optional[int] = 2,
    max_iter_ACTION: Optional[int] = 50,
    min_specificity_z_threshold: Optional[int] = -3,
    distance_metric: Optional[str] = "jsd",
    nn_approach: Optional[str] = "k*nn",
    network_density: Optional[int] = 1,
    network_k: Optional[int] = 30,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 1000,
    layout_algorithm: Optional[str] = "tumap",
    layout_in_parallel: Optional[bool] = True,
    unification_violation_threshold: Optional[float] = 0,
    footprint_key: Optional[str] = "archetype_footprint",
    footprint_alpha: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False
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
    else:
        S_r = adata.obsm[reduction_key].astype(dtype=np.float64)

    mr = decomp.ACTIONMR(
        k_max,
        n_iter=max_iter_ACTION,
        tol=1e-300,
        min_specificity_z_threshold=min_specificity_z_threshold,
        min_cells_per_arch=min_cells_per_archetype,
        unification_violation_threshold=unification_violation_threshold,
        thread_no=thread_no,
    )

    H_unified = mr.fit_transform(S_r)
    C_unified = mr.coeff_.T

    adata.obsm["C_unified"] = C_unified
    adata.obsm["H_unified"] = H_unified
    adata.uns.setdefault("obsm_annot", {}).update(
        {
            "C_unified": {"type": np.array([b"internal"], dtype=object)},
            "H_unified": {"type": np.array([b"internal"], dtype=object)},
        }
    )

    adata.obs["assigned_archetype"] = pd.Categorical(
        values=mr.assigned_archetype.astype(int),
        categories=natsorted(map(int, np.unique(mr.assigned_archetype))),
    )

    C_stacked = mr.stacked_coeffs
    H_stacked = mr.stacked_loadings
    adata.obsm["C_stacked"] = C_stacked.T
    adata.obsm["H_stacked"] = H_stacked
    adata.uns.setdefault("obsm_annot", {}).update(
        {
            "C_stacked": {"type": np.array([b"internal"], dtype=object)},
            "H_stacked": {"type": np.array([b"internal"], dtype=object)},
        }
    )

    # Build ACTIONet
    net.build(
        data=adata,
        algorithm=nn_approach,
        distance_metric=distance_metric,
        density=network_density,
        mutual_edges_only=mutual_edges_only,
        k=network_k,
        net_key=net_key,
        data_key="H_stacked",
        thread_no=thread_no,
        copy=False,
        return_raw=False
    )

    # Uses `reduction_key` for initialization of coordinates
    net.layout(
        data=adata,
        algorithm=layout_algorithm,
        initial_coordinates=None,
        compactness_level=layout_compactness,
        n_epochs=layout_epochs,
        thread_no=thread_no if layout_in_parallel else 1,
        reduction_key=reduction_key,
        net_key=net_key,
        seed=seed,
        copy=False,
        return_raw=False
    )

    # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
    net.compute_core_centrality(
        adata=adata,
        G=None,
        assignments=None,
        net_key=net_key,
        assignment_key="assigned_archetype",
        copy=False,
        return_raw=False
    )

    # Smooth archetype footprints
    net.diffusion(
        adata=adata,
        G=None,
        H_unified=None,
        net_key=net_key,
        footprint_key=footprint_key,
        alpha_val=footprint_alpha,
        thread_no=thread_no,
        copy=False,
        return_raw=False
    )

    # Compute gene specificity for each archetype
    pp.compute_archetype_feature_specificity(
        adata=adata,
        S=S,
        H=None,
        layer_key=None,
        footprint_key=footprint_key,
        copy=False,
        return_raw=False
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
