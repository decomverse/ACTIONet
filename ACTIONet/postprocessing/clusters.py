from typing import Optional, Union, Literal, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from .. import tools as tl


def __compute_feature_specificity(S, sample_assignments, thread_no=0):
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        out = _an.compute_cluster_feature_specificity(S, sample_assignments, thread_no)
    else:
        S = np.array(S, dtype=np.float64)
        out = _an.compute_cluster_feature_specificity_full(
            S, sample_assignments, thread_no
        )

    return out


def feature_specificity(
    data: Union[np.ndarray, sparse.spmatrix] = None,
    cluster_attr: Union[str, list, pd.Series, np.ndarray, None] = None,
    output_prefix: Optional[str] = None,
    layer_key: Optional[str] = None,
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of clusters

    Uses cluster membership vector to estimate markers (disjoint clustering)
    Parameters
    ----------
    data
        AnnData object or expression matrix of shape `n_obs` × `n_vars`..
    cluster_attr
        List-like object of length 'data.shape[0]' containing clusters labels or key in `adata.obs` that holds the clustering variable.
    output_prefix
        String to prefix keys in 'adata.varm' where output is stored.
    layer_key
        Key of 'layers' to containing gene expression matrix. Defaults to 'adata.X' if 'None'.
    thread_no : Optional[int], optional
        Number of threads.
    copy
        Return a copy instead of writing to adata.
    return_raw
        If 'adata' is given, return dict of 'compute_cluster_feature_specificity()' output instead of storing to 'adata'.

    Returns
    -------
        adata : anndata.AnnData
        If `copy=True` returns None or else adds fields to `adata`:
        `.varm[f'{cluster_key}_feature_specificity']`

        specificity: dict
        If 'return_raw=True', returns dict containing 'average_profile', 'upper_significance', and 'lower_significance' matrices.
    """

    data_is_AnnData = isinstance(data, AnnData)

    if cluster_attr is None:
        raise ValueError(f"'cluster_attr' cannot be 'None.")

    if data_is_AnnData:
        adata = data.copy() if copy else data

        if layer_key is not None:
            if layer_key not in adata.layers.keys():
                raise ValueError(f"'{layer_key}' not in 'adata.layers.keys()'.")
            S = adata.layers[layer_key]
        else:
            S = adata.X

        sa = tl.get_data_or_split(adata=adata, attr=cluster_attr, to_return="levels")
        clusters = sa["index"]

    else:
        adata = None
        if len(cluster_attr) != data.shape[0]:
            raise ValueError(f"'len(cluster_attr)' must equal 'data.shape[0]'.")
        S = data
        clusters = pd.factorize(list(cluster_attr), sort=True)[0]

    S = S.T.astype(dtype=np.float64)

    specificity_out = __compute_feature_specificity(
        S=S, sample_assignments=clusters + 1, thread_no=thread_no
    )

    if return_raw or not data_is_AnnData:
        return specificity_out
    else:
        adata.varm[f"{output_prefix}_feature_specificity"] = specificity_out[
            "upper_significance"
        ]
        adata.uns.setdefault("varm_annot", {}).update(
            {
                f"{output_prefix}_feature_specificity": {
                    "type": np.array([b"reduction"], dtype=object)
                },
            }
        )

        return adata if copy else None


def annotate(
    adata: AnnData,
    markers=None,
    labels=None,
    scores=None,
    thread_no: Optional[int] = 0,
    cluster_key: Optional[str] = "leiden",
    specificity_key: Optional[str] = "leiden_feature_specificity",
):

    if markers is not None:
        feature_names = pd.Series(
            [
                x.decode() if isinstance(x, (bytes, bytearray)) else x
                for x in list(adata.var.index)
            ]
        )

        marker_mat = sparse.csc_matrix(
            pd.DataFrame(
                [feature_names.isin(markers[key]) * 1 for key in markers.keys()]
            ).T
        )
        spec_mat = adata.varm[specificity_key]

        assessment_out = _an.assess_enrichment(
            scores=spec_mat, associations=marker_mat, thread_no=thread_no
        )

        logPvals = assessment_out["logPvals"].T
        Enrichment = pd.DataFrame(
            logPvals, columns=markers.keys(), index=np.arange(1, logPvals.shape[0] + 1),
        )
        annotations = pd.Series(markers.keys())
        idx = np.argmax(logPvals, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(logPvals, axis=1)

    elif labels is not None:
        cluster_dict = tl.get_data_or_split(
            adata=adata, attr=cluster_key, to_return="levels"
        )
        X1 = np.array(
            pd.DataFrame(
                [
                    (cluster_dict["index"] == k) * 1
                    for k in range(len(cluster_dict["keys"]))
                ]
            ).T
        )

        labels_dict = tl.get_data_or_split(adata=adata, attr=labels, to_return="levels")
        X2 = np.array(
            pd.DataFrame(
                [
                    (labels_dict["index"] == k) * 1
                    for k in range(len(labels_dict["keys"]))
                ]
            ).T
        )

        XI = _an.XICOR(X1, X2)
        Z = np.sign(tl.scale_matrix(X1).T @ tl.scale_matrix(X2)) * XI["Z"]
        Enrichment = pd.DataFrame(
            Z, index=cluster_dict["keys"], columns=labels_dict["keys"]
        )

        annotations = pd.Series(labels_dict["keys"])
        idx = np.argmax(Z, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(Z, axis=1)

    elif scores is not None:
        cluster_dict = tl.get_data_or_split(
            adata=adata, attr=cluster_key, to_return="levels"
        )
        X1 = np.array(
            pd.DataFrame(
                [
                    (cluster_dict["index"] == k) * 1
                    for k in range(len(cluster_dict["keys"]))
                ]
            ).T
        )

        if isinstance(scores, str) or len(scores) == 1:
            if scores in adata.obs.keys():
                scores = pd.DataFrame(adata.obs[scores])
            elif scores in adata.obsm.keys():
                scores = pd.DataFrame(adata.obsm[scores])
            else:
                raise ValueError("label attribute %s not found" % scores)

        X2 = np.array(scores)

        XI = _an.XICOR(X1, X2)
        Z = np.sign(tl.scale_matrix(X1).T @ tl.scale_matrix(X2)) * XI["Z"]
        Enrichment = pd.DataFrame(Z, index=cluster_dict["keys"], columns=scores.columns)

        annotations = scores.columns
        idx = np.argmax(Z, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(Z, axis=1)

    return Label, Confidence, Enrichment

