from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def centrality(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: str = "local_coreness",
    annotations: Union[np.ndarray, list, pd.Series] = None,
    annotations_key: Optional[str] = "assigned_archetype",
    net_key: Optional[str] = "ACTIONet",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes node centrality scores

    Compute node centralities using different measures

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    algorithm: str
        centrality algorithm. Can be "coreness", "local_coreness", "pagerank", "personalized_pagerank", default is "local_coreness"
    annotations:
        list-like object containing sample annotations/scores of each observation (for localized measures).
        Required if 'adata=None'.
    annotations_key:
        Key of 'adata.obs' containing list-like object of sample annotations/scores (default="assigned_archetype").
        Ignored if data is not an AnnData object.
    net_key:
        Key of 'adata.obsp' containing adjacency matrix (default="ACTIONet").
        Ignored if data is not an AnnData object.
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
    if data is not None:
        if isinstance(data, AnnData):
            adata = data.copy() if copy else data
            annotations = (
                annotations if annotations is not None else adata.obsm[annotations_key]
            )
            if net_key in adata.obsp.keys():
                G = adata.obsp[net_key]
            else:
                raise Exception("missing %s in adata.obsp of AnnData" % net_key)
        else:
            G = data

    if G is None or annotations is None:
        raise ValueError("'G' and 'assignments' cannot be NoneType.")
    if not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
    if not isinstance(annotations, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'assignments' must be numpy.ndarray or sparse.spmatrix.")

    if algorithm == "local_coreness":
        node_centrality = _an.compute_archetype_core_centrality(G, annotations)
    elif algorithm == "coreness":
        node_centrality = _an.compute_core_number(G)
    # elif algorithm == "pagerank":
    #     u=np.empty(G.shape[0]);
    #     u.fill(1/G.shape[0])
    #     pr = diffusion(G, u)
    #     node_centrality = pr[:, 1]

    if return_raw or adata is None:
        return node_centrality
    else:
        adata.obs["node_centrality"] = node_centrality
        return adata if copy else None


def diffusion(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    scores: Union[np.ndarray, sparse.spmatrix] = None,
    algorithm: Optional[str] = "pagerank_chebyshev",
    alpha_val: Optional[float] = 0.85,
    max_it: Optional[int] = 5,
    threshold: Optional[float] = 1e-8,
    thread_no: Optional[int] = 0,
    net_key: Optional[str] = "ACTIONet",
    scores_key: Optional[str] = "H_stacked",
    smoothed_scores_key: Optional[str] = "archetype_footprint",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes smoothed scores using network diffusion

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    scores : Union[np.ndarray, sparse.spmatrix], optional
        Input scores, by default None
    algorithm : Optional[str], optional
        Diffusion algorithm to use. Can be "pagerank", "pagerank_sym", by default "pagerank_chebyshev"
    alpha_val : Optional[float], optional
        Diffusion parameter. Larger values results in more long-range diffusion, by default 0.85
    max_it : Optional[int], optional
        [description], by default 5
    threshold : Optional[float], optional
        [description], by default 1e-8
    thread_no : Optional[int], optional
        Number of threads to use, by default 0
    net_key : Optional[str], optional
        Key of 'adata.obsp' containing adjacency matrix to use. (default="ACTIONet")
        Ignored if 'adata=None'.
    scores_key : Optional[str], optional
        Key of 'adata.obsm' containing scores. (default="H_stacked")
        Ignored if `adata=None`
    smoothed_scores_key : Optional[str], optional
        Key of 'adata.obsm' to store smoothed scores. (default="archetype_footprint")
        Ignored if `adata=None`
    copy : Optional[bool], optional
        If 'adata' is given, return a copy instead of writing to `adata` (default=False)
    return_raw : Optional[bool], optional
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata' (default=False)

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm["archetype_footprint"]`

    smoothed_scores : np.ndarray
        If `adata=None` or `return_raw=True`, returns array of archetype footprint.
    """
    if data is not None:
        if isinstance(data, AnnData):
            adata = data.copy() if copy else data
            scores = scores if scores is not None else adata.obsm[scores_key]
            if net_key in adata.obsp.keys():
                G = adata.obsp[net_key]
            else:
                raise Exception("missing %s in adata.obsp of AnnData" % net_key)
        else:
            G = data

    if G is None or scores is None:
        raise ValueError("'G' and 'scores' cannot be NoneType.")
    if not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
    if not isinstance(scores, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'scores' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    scores = scores.astype(dtype=np.float64)

    if not sparse.issparse(scores):
        scores = sparse.csc_matrix(scores)

    smoothed_scores = _an.compute_network_diffusion_fast(
        G=G, X0=scores, alpha=alpha_val, thread_no=thread_no, max_it=max_it
    )

    smoothed_scores = np.array(smoothed_scores, dtype=np.float64)

    if return_raw or adata is None:
        return smoothed_scores
    else:
        adata.obsm[smoothed_scores_key] = smoothed_scores
        return adata if copy else None


def construct_backbone(
    adata: Optional[AnnData] = None,
    footprint: Union[np.ndarray, sparse.spmatrix] = None,
    net_key: Optional[str] = "ACTIONet",
    footprint_key: Optional[str] = "archetype_footprint",
    scale: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 100,
    alpha_val: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    max_it: Optional[int] = 5,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
        else:
            raise ValueError("'adata' is not an AnnData object.")

    if footprint_key not in adata.obsm.keys():
        diffusion(
            adata=adata,
            G=None,
            H_unified=None,
            net_key=net_key,
            footprint_key=footprint_key,
            alpha_val=alpha_val,
            thread_no=thread_no,
            max_it=max_it,
            copy=False,
            return_raw=False,
        )

    footprint = footprint if footprint is not None else adata.obsm[footprint_key]

    if scale:
        W = np.exp(
            (footprint - np.mean(footprint, axis=0)) / np.std(footprint, axis=0, ddof=1)
        )
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
        seed=seed,
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
        {"backbone": backbone,}
    )

    return adata if copy else None


def correct_cell_annotations(
    adata: Optional[AnnData],
    initial_labels: Optional[list],
    iters: Optional[int] = 3,
    lambda_param: Optional[float] = 0,
    sig_threshold: Optional[int] = 3,
    min_cell_fraction: Optional[float] = 0.001,
):
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

