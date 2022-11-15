from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
import ACTIONet as an
from ACTIONet.network.cluster import cluster as anet_cluster
from ACTIONet.network.propagation import propagate
from ACTIONet.postprocessing.clusters import feature_specificity
from scipy.sparse import csc_matrix


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
    alpha: Optional[float] = 0.85,
    post_correction: Optional[bool] = False,
    thread_no: Optional[int] = 0,
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
    S = sparse.csc_matrix(adata.X.T)
    marker_mat = sparse.csc_matrix(masks)

    G = sparse.csc_matrix(adata.obsp[network_key])

    algorithm = "parametric"  # For now!

    if algorithm == "parametric":
        out = _an.aggregate_genesets(
            G, S, marker_mat, network_normalization_code, alpha, thread_no
        )
        marker_stats = out["stats_norm_smoothed"]
        Enrichment = pd.DataFrame(
            marker_stats, index=adata.obs.index, columns=markers.keys()
        )

        marker_stats_raw = out["stats_norm"]
        marker_stats_raw[marker_stats_raw < 0] = 0
        G_norm = csc_matrix(_an.normalize_spmat(G, 1).T)
        logPvals = _an.assess_label_enrichment(G_norm, marker_stats_raw)

        annotations = pd.Series(markers.keys())
        idx = np.argmax(logPvals, axis=1)
        Label = annotations[idx]
        Label.index = adata.obs.index
        Confidence = np.max(logPvals, axis=1)

    if post_correction == True:
        Label = correct_labels(adata=adata, initial_labels=Label, return_raw=True)

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

    if not isinstance(adata, AnnData):
        raise Exception("adata must be an AnnData object")

    adata = adata.copy() if copy else adata

    anet_cluster(
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

    feature_specificity(
        data=adata,
        cluster_attr=final_clusters_key,
        output_prefix=final_clusters_key,
        thread_no=thread_no,
    )

    return adata if copy else None


def infer_missing_labels(
    adata: AnnData,
    initial_labels: Union[str, np.ndarray, list, pd.Series],
    algorithm: str = "lpa",
    lambda_param: int = 0,
    iters: int = 3,
    sig_threshold: float = 3,
    thread_no: int = 0,
    output_key: str = "inferred_labels",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> pd.Series:

    if not isinstance(adata, AnnData):
        raise Exception("adata must be an AnnData object")

    adata = adata.copy() if copy else adata

    if isinstance(initial_labels, str):
        if initial_labels in adata.obs.keys():
            labels = pd.Series(
                adata.obs[initial_labels].astype("str"),
                index=adata.obs.index.values.astype("str"),
            )
        else:
            raise ValueError("labels attribute %s not found" % initial_labels)
    else:
        if isinstance(initial_labels, list):
            labels = pd.Series(
                [str(x) for x in initial_labels],
                index=adata.obs.index.values.astype("str"),
            )
        elif isinstance(initial_labels, pd.Series):
            labels = pd.Series(
                initial_labels.astype("str"), index=adata.obs.index.values.astype("str")
            )

    fixed_samples = np.where(labels != "nan")[0]

    updated_labels = propagate(
        data=adata,
        labels=labels,
        fixed_samples=fixed_samples,
        return_raw=True,
        algorithm=algorithm,
        lambda_param=lambda_param,
        iters=iters,
        sig_threshold=sig_threshold,
        thread_no=thread_no,
    )

    if return_raw or not isinstance(adata, AnnData):
        return updated_labels
    else:
        adata.obsm[output_key] = updated_labels
        return adata if copy else None


def correct_labels(
    adata: AnnData,
    initial_labels: Union[str, np.ndarray, list, pd.Series],
    algorithm: str = "lpa",
    lambda_param: int = 0,
    iters: int = 3,
    sig_threshold: float = 3,
    thread_no: int = 0,
    output_key: str = "inferred_labels",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> pd.Series:

    if not isinstance(adata, AnnData):
        raise Exception("adata must be an AnnData object")

    if initial_labels is None:
        raise Exception("initial_labels must be provided")

    adata = adata.copy() if copy else adata

    if isinstance(initial_labels, str):
        if initial_labels in adata.obs.keys():
            labels = pd.Series(
                adata.obs[initial_labels].astype("str"),
                index=adata.obs.index.values.astype("str"),
            )
        else:
            raise ValueError("labels attribute %s not found" % initial_labels)
    else:
        if isinstance(initial_labels, list):
            labels = pd.Series(
                [str(x) for x in initial_labels],
                index=adata.obs.index.values.astype("str"),
            )
        elif isinstance(initial_labels, pd.Series):
            labels = pd.Series(
                initial_labels.astype("str"), index=adata.obs.index.values.astype("str")
            )
        elif isinstance(initial_labels, np.ndarray):
            labels = pd.Series(
                initial_labels.astype("str"), index=adata.obs.index.values.astype("str")
            )

    updated_labels = an.nt.propagate(
        data=adata,
        labels=labels,
        fixed_samples=[],
        return_raw=True,
        algorithm=algorithm,
        lambda_param=lambda_param,
        iters=iters,
        sig_threshold=sig_threshold,
        thread_no=thread_no,
    )

    if return_raw or not isinstance(adata, AnnData):
        return updated_labels
    else:
        adata.obs[output_key] = updated_labels
        return adata if copy else None
