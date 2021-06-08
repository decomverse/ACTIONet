from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def compute_archetype_core_centrality(
    adata: Optional[AnnData] = None,
    G: Union[np.ndarray, sparse.spmatrix] = None,
    assignments: Union[np.ndarray, list, pd.Series] = None,
    net_name: Optional[str] = "ACTIONet",
    assignment_name: Optional[str] = "assigned_archetype",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False
) -> Union[AnnData, np.ndarray, None]:
    """
    Computes node centrality scores

    Uses graph core-ness to compute node centralities

    Parameters
    ----------
    adata
        AnnData object possibly containing 'assignment_name' in '.obs' and 'net_name' in '.obsp'.
    G:
        Adjacency matrix to use for computing centrality.
        Required if 'adata=None'.
    assignments:
        list-like object containing archetype assignments of each observation.
        Required if 'adata=None'.
    assignment_name:
        Key of 'adata.obs' containing list-like object of archetype assignments (default="assigned_archetype").
        Ignored if 'adata=None'.
    net_name:
        Key of 'adata.obsp' containing adjacency matrix to use for 'G' in 'compute_archetype_core_centrality()' (default="ACTIONet").
        Ignored if 'adata=None'.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata'.

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obs['node_centrality']`

    node_centrality : np.ndarray
        If 'adata=None' or 'return_raw=True', returns array of node centrality scores for each observation.
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            G = G if G is not None else adata.obsp[net_name]
            assignments = assignments if assignments is not None else adata.obs[assignment_name]
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if G is None or assignments is None:
            raise ValueError("'G' and 'S_r' cannot be NoneType if 'adata=None'.")
        if not isinstance(G, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(assignments, (np.ndarray, list, pd.Series)):
            raise ValueError("'S_r' must be list, numpy.ndarray or pandas.Series.")

    node_centrality = _an.compute_archetype_core_centrality(G, assignments)

    if return_raw or adata is None:
        return node_centrality
    else:
        adata.obs["node_centrality"] = node_centrality
        return adata if copy else None


def compute_network_diffusion(
    adata: Optional[AnnData] = None,
    G: Union[np.ndarray, sparse.spmatrix] = None,
    H_unified: Union[np.ndarray, sparse.spmatrix] = None,
    net_name: Optional[str] = "ACTIONet",
    footprint_key: Optional[str] = "archetype_footprint",
    alpha_val: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False
) -> Union[AnnData, np.ndarray, None]:

    """
    Computes archetype footprint via network diffusion.

    Parameters
    ----------
    adata
        AnnData object possibly containing '.obsm["H_unified]' and 'net_name' in '.obsp'.
    G:
        Adjacency matrix to use for computing diffusion.
        Required if 'adata=None'.
    H_unified:
        Matrix containing output 'H_unified' of 'unify_archetypes()' to use for computing diffusion.
        Required if 'adata=None', otherwise retrieved from '.obsm["H_unified"]'
    net_name:
        Key of 'adata.obsp' containing adjacency matrix to use for 'G' in 'compute_archetype_core_centrality()' (default="ACTIONet").
        Ignored if 'adata=None'.
    footprint_key:
        Key of 'adata.obsp' to store archetype footprint (default="archetype_footprint").
        Ignored if 'adata=None'.
    alpha_val:
        Diffusion parameter between 0-1.
    thread_no:
        Number of threads. Defaults to number of threads available - 2.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata'.

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm["archetype_footprint"]`

    archetype_footprint : np.ndarray
        If 'adata=None' or 'return_raw=True', returns array of archetype footprint.
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            G = G if G is not None else adata.obsp[net_name]
            H_unified = H_unified if H_unified is not None else adata.obsm["H_unified"]
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if G is None or H_unified is None:
            raise ValueError("'G' and 'H_unified' cannot be NoneType if 'adata=None'.")
        if not isinstance(G, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(H_unified, (np.ndarray, list, pd.Series)):
            raise ValueError("'H_unified' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    H_unified = H_unified.astype(dtype=np.float64)

    if not sparse.issparse(H_unified):
        H_unified = sparse.csc_matrix(H_unified)

    archetype_footprint = _an.compute_network_diffusion(G, H_unified, alpha=alpha_val, thread_no=thread_no)
    archetype_footprint = np.array(archetype_footprint, dtype=np.float64)

    if return_raw or adata is None:
        return archetype_footprint
    else:
        adata.obsm[footprint_key] = archetype_footprint
        return adata if copy else None


def construct_backbone(
    adata: Optional[AnnData] = None,
    footprint: Union[np.ndarray, sparse.spmatrix] = None,
    net_name: Optional[str] = "ACTIONet",
    footprint_key: Optional[str] = "archetype_footprint",
    scale: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 100,
    alpha_val: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
        else:
            raise ValueError("'adata' is not an AnnData object.")

    if footprint_key not in adata.obsm.keys():
        compute_network_diffusion(
            adata=adata,
            G=None,
            H_unified=None,
            net_name=net_name,
            footprint_key=footprint_key,
            alpha_val=alpha_val,
            thread_no=thread_no,
            copy=False,
            return_raw=False
        )

    footprint = footprint if footprint is not None else adata.obsm[footprint_key]

    if scale:
        W = np.exp( (footprint - np.mean(footprint, axis=0)) / np.std(footprint, axis=0, ddof=1) )
    else:
        W = np.exp(footprint)

    W = sparse.csc_matrix(W)
    
    arch_vis_out = _an.transform_layout(
        W=W,
        coor2D=adata.obsm["ACTIONet2D"].T,
        coor3D=adata.obsm["ACTIONet3D"].T,
        colRGB=adata.obsm["denovo_color"].T,
        compactness_level=layout_compactness,
        n_epochs=layout_epochs,
        thread_no=thread_no,
        seed=seed
    )
    
    arch_G = _an.compute_full_sim(footprint, thread_no)
    np.fill_diagonal(arch_G, 0)

    backbone = {
        "G": arch_G,
        "coordinates": arch_vis_out["coordinates"],
        "coordinates_3D": arch_vis_out["coordinates_3D"],
    }

    # adata.uns['metadata']["ACTIONet_backbone"] = backbone
    adata.uns.setdefault("metadata", {}).update(
        {
            "backbone": backbone,
        })

    return adata if copy else None


def correct_cell_annotations(adata: Optional[AnnData],
                             initial_labels: Optional[list],
                             iters: Optional[int] = 3,
                             lambda_param: Optional[float] = 0,
                             sig_threshold: Optional[int] = 3,
                             min_cell_fraction: Optional [float] = 0.001):
    """
    Uses a variant of the label propagation algorithm to correct likely noisy labels

    #' @param ace Input results to be clustered
    #' (alternatively it can be the ACTIONet igraph object)
    #' @param initial_labels Annotations to correct with missing values (NA) in it.
    #' It can be either a named annotation (inside ace$annotations) or a label vector.
    #' @param LFR.threshold How aggressively to update labels. The smaller the value, the more labels will be changed (default=2)
    #' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
    #' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
    #' @param min.cell.fraction Annotations with less that this fraction will be removed
    #'
    #' @return ace with updated annotations added to ace$annotations
    #'
    #' @examples
    #' ace = add.cell.annotations(ace, cell.labels, 'input_annotations')
    #' ace = correct.cell.annotations(ace, 'input_annotations', 'updated_annotations')
    #' @export
    """
    pass

def enhance_cell_annotations():
    pass

