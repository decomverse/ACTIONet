from pickle import TRUE
from typing import Optional, Union
from typing_extensions import Literal

import numpy as np
import pandas as pd
from scipy import sparse
from anndata import AnnData

import _ACTIONet as _an
import ACTIONet as an


def filter(
    adata: AnnData,
    centrality_key: Optional[str] = "node_centrality",
    z_threshold: Optional[float] = -1.65,
    output_key: Optional[str] = "is_filtered",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
):
    z = _an.normalize_mat(adata.obs[centrality_key], -1)
    is_filtered = z < z_threshold

    if return_raw:
        return is_filtered
    else:
        adata.obs[output_key] = is_filtered

    return adata if copy else None


def annotate(
    adata: AnnData,
    markers,
    algorithm: Optional[str] = "parametric",
    network_normalization_method: Optional[str] = "pagerank_sym",
    network_key: Optional[str] = "ACTIONet",
    prenorm: Optional[str] = 1,
    gene_scaling_method: Optional[int] = 0,
    pre_alpha: Optional[float] = 0.15,
    post_alpha: Optional[float] = 0.9,
    perm_no: Optional[int] = 30,
):
    network_normalization_code = 0
    if network_normalization_method == "pagerank_sym":
        network_normalization_code = 2

    feature_names = pd.Series(
        [
            x.decode() if isinstance(x, (bytes, bytearray)) else x
            for x in list(adata.var.index)
        ]
    )
    masks = np.array(
        pd.DataFrame([feature_names.isin(markers[key]) * 1 for key in markers.keys()]).T
    )
    selected_features_mask = np.sum(masks, axis=1) > 0
    sub_S = sparse.csc_matrix(adata.X[:, selected_features_mask].T)
    marker_mat = sparse.csc_matrix(masks[selected_features_mask, :])

    G = sparse.csc_matrix(adata.obsp[network_key])

    if algorithm == "parametric":
        marker_stats = _an.aggregate_genesets_weighted_enrichment_permutation(
            G,
            sub_S,
            marker_mat,
            network_normalization_method=network_normalization_code,
            expression_normalization_method=prenorm,
            gene_scaling_method=gene_scaling_method,
            pre_alpha=pre_alpha,
            post_alpha=post_alpha,
            perm_no=perm_no,
        )
    elif algorithm == "nonparametric":
        marker_stats = _an.aggregate_genesets_weighted_enrichment(
            G,
            sub_S,
            marker_mat,
            network_normalization_method=network_normalization_code,
            expression_normalization_method=prenorm,
            gene_scaling_method=gene_scaling_method,
            pre_alpha=pre_alpha,
            post_alpha=post_alpha,
        )

    Enrichment = pd.DataFrame(
        marker_stats, index=adata.obs.index, columns=markers.keys()
    )
    annotations = pd.Series(markers.keys())
    idx = np.argmax(marker_stats, axis=1)
    Label = annotations[idx]
    Label.index = adata.obs.index
    Confidence = np.max(marker_stats, axis=1)
    return Label, Confidence, Enrichment


def cluster(
    adata: AnnData,
    algorithm: str = "leiden",
    resolution: Optional[float] = 1.0,
    initial_clusters: Optional[Union[np.ndarray, list, pd.Series]] = None,
    initial_clusters_key: Optional[str] = "assigned_archetype",
    final_clusters_key: Optional[str] = "leiden",
    seed: Optional[int] = 0,
    net_key: Optional[str] = "ACTIONet",
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Cluster cells

    Cluster cells

    Parameters
    ----------
    adata : AnnData
        Adjacency matrix of the input graph or AnnData object containing the network.
    algorithm: str
        centrality algorithm. Can be "Leiden", "fix", default is "Leiden"
    resolution: float
        Resolution of the Leiden clustering. Larger values results in more clusters.
    initial_clusters:
        Initial clusters.
    initial_clusters_key:
        Key of 'adata.obs' containing the initial clusters (default="assigned_archetype").
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

    data_is_AnnData = isinstance(adata, AnnData)
    if not data_is_AnnData:
        raise Exception("adata must be an AnnData object")

    adata = adata.copy() if copy else adata

    an.net.cluster(
        data=adata,
        algorithm=algorithm,
        resolution=resolution,
        initial_clusters=initial_clusters,
        initial_clusters_key=initial_clusters_key,
        final_clusters_key=final_clusters_key,
        seed=seed,
        net_key=net_key,
        copy=False,
    )

    an.po.clusters.feature_specificity(
        data=adata,
        cluster_attr=final_clusters_key,
        output_prefix=final_clusters_key,
        thread_no=thread_no,
    )

    return adata if copy else None


def infer_missing_labels(
    adata: AnnData,
    algorithm: str = "lpa",
    initial_labels: Union[str, np.ndarray, list, pd.Series] = None,
    lambda_param: int = 0,
    iters: int = 3,
    sig_threshold: float = 3,
) -> pd.Series:

    data_is_AnnData = isinstance(adata, AnnData)
    if not data_is_AnnData:
        raise Exception("adata must be an AnnData object")

    if isinstance(initial_labels, str):
        if initial_labels in adata.obs.keys():
            labels = pd.DataFrame(adata.obs[initial_labels])
        else:
            raise ValueError("labels attribute %s not found" % initial_labels)
    else:
        labels = pd.Series(initial_labels)

    fixed_samples = np.where(labels != None)

    updated_labels = an.net.propagate(
        data=adata,
        labels=labels,
        fixed_samples=fixed_samples,
        return_raw=True,
        algorithm=algorithm,
        lambda_param=lambda_param,
        iters=iters,
        sig_threshold=sig_threshold,
    )

    return updated_labels


def correct_labels(
    adata: AnnData,
    algorithm: str = "lpa",
    initial_labels: Union[str, np.ndarray, list, pd.Series] = None,
    lambda_param: int = 0,
    iters: int = 3,
    sig_threshold: float = 3,
) -> pd.Series:

    data_is_AnnData = isinstance(adata, AnnData)
    if not data_is_AnnData:
        raise Exception("adata must be an AnnData object")

    if isinstance(initial_labels, str):
        if initial_labels in adata.obs.keys():
            labels = pd.DataFrame(adata.obs[initial_labels])
        else:
            raise ValueError("labels attribute %s not found" % initial_labels)
    else:
        labels = pd.Series(initial_labels)

    updated_labels = an.net.propagate(
        data=adata,
        labels=labels,
        fixed_samples=None,
        return_raw=True,
        algorithm=algorithm,
        lambda_param=lambda_param,
        iters=iters,
        sig_threshold=sig_threshold,
    )

    return updated_labels
