from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from ACTIONet.network.build import build
from ACTIONet.tools.utils_public import double_normalize, scale_matrix


def __compute_feature_specificity(S, H, thread_no=0):
    H = np.array(H, dtype=np.float64)
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        out = _an.compute_archetype_feature_specificity(S, H, thread_no)
    else:
        S = np.array(S, dtype=np.float64)
        out = _an.compute_archetype_feature_specificity_full(S, H, thread_no)

    return out


def feature_specificity(
    adata: Optional[AnnData] = None,
    S: Union[np.ndarray, sparse.spmatrix] = None,
    H: Union[np.ndarray, sparse.spmatrix] = None,
    layer_key: Optional[str] = None,
    footprint_key: Optional[str] = "archetype_footprint",
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, dict, None]:
    """Computes Feature (i.e., gene) specificity of archetypes \
    Uses Archetype footprints to estimate markers (soft clustering)
    :param adata: AnnData object possibly containing '.layers["layer_key"]' and '.obsm["footprint_key"]'.
    :param S: `n_obs` × `n_vars` gene expression matrix. \
        Required if 'adata=None', otherwise retrieved from '.layers["layer_key"]' or '.X' if 'layer_key=None'.
    :param H: Matrix-like object containing archetype footprint \
        Required if 'adata=None', otherwise retrieved from '.obsm["footprint_key"]'
    :param layer_key: Key of 'layers' to retrieve gene expression matrix.
    :param footprint_key: Key in `adata.obsm` that holds the archetype footprint.
    :param thread_no: Number of threads.
    :param copy: If 'adata' is given, return a copy instead of writing to `adata`
    :param return_raw: If 'adata' is given, return dict of 'compute_archetype_feature_specificity()' output instead of storing to 'adata'.
    ...
    :return adata : anndata.AnnData \
        If `copy=True` returns None or else adds fields to `adata` \
        `.varm["unified_feature_profile"]` \
        `.varm["unified_feature_specificity"]`
    :return specificity: dict \
        If 'adata=None' or 'return_raw=True', returns dict containing 'unified_feature_profile' and 'unified_feature_specificity' matrices.
    """
    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            if layer_key is not None:
                S = S.T if S is not None else adata.layers[layer_key].T
            else:
                S = S.T if S is not None else adata.X.T
            H = H.T if H is not None else adata.obsm[footprint_key].T
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if S is None or H is None:
            raise ValueError("'S' and 'H' cannot be NoneType if 'adata=None'.")
        if not isinstance(S, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'S' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(H, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'H' must be numpy.ndarray or sparse.spmatrix.")
        S = S.T
        H = H.T

    S = S.astype(dtype=np.float64)
    H = H.astype(dtype=np.float64)

    specificity = __compute_feature_specificity(S=S, H=H, thread_no=thread_no)

    if return_raw or adata is None:
        return specificity
    else:
        adata.varm["unified_feature_profile"] = specificity["archetypes"]
        adata.varm["unified_feature_specificity"] = specificity["upper_significance"]

        adata.uns.setdefault("varm_annot", {}).update(
            {
                "unified_feature_profile": {"type": np.array([b"internal"], dtype=object)},
                "unified_feature_specificity": {"type": np.array([b"reduction"], dtype=object)},
            }
        )

        return adata if copy else None


def annotate(
    adata: AnnData,
    markers=None,
    labels=None,
    scores=None,
    thread_no: Optional[int] = 0,
    archetype_key: Optional[str] = "H_unified",
    specificity_key: Optional[str] = "unified_feature_specificity",
):

    if markers is not None:
        feature_names = pd.Series([x.decode() if isinstance(x, (bytes, bytearray)) else x for x in list(adata.var.index)])

        marker_mat = sparse.csc_matrix(pd.DataFrame([feature_names.isin(markers[key]) * 1 for key in markers.keys()]).T)
        spec_mat = adata.varm[specificity_key]

        assessment_out = _an.assess_enrichment(scores=spec_mat, associations=marker_mat, thread_no=thread_no)

        logPvals = assessment_out["logPvals"].T
        Enrichment = pd.DataFrame(
            logPvals,
            columns=markers.keys(),
            index=np.arange(1, logPvals.shape[0] + 1),
        )
        annotations = pd.Series(markers.keys())
        idx = np.argmax(logPvals, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(logPvals, axis=1)

    elif labels is not None:
        X1 = adata.obsm["H_unified"].toarray()

        if isinstance(adata, AnnData) and isinstance(labels, str):
            labels = adata.obs[labels]

        labels_index, labels_keys = pd.factorize(labels, sort=True)
        X2 = np.array(pd.DataFrame([(labels_index == k) * 1 for k in range(len(labels_keys))]).T)

        XI = _an.XICOR(X1, X2)
        Z = np.sign(scale_matrix(X1).T @ scale_matrix(X2)) * XI["Z"]
        Enrichment = pd.DataFrame(Z, index=np.arange(1, X1.shape[1] + 1), columns=labels_keys)

        annotations = pd.Series(labels_keys)
        idx = np.argmax(Z, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(Z, axis=1)

    elif scores is not None:
        X1 = adata.obsm["H_unified"].toarray()

        if isinstance(scores, str) or len(scores) == 1:
            if scores in adata.obs.keys():
                scores = pd.DataFrame(adata.obs[scores])
            elif scores in adata.obsm.keys():
                scores = pd.DataFrame(adata.obsm[scores])
            else:
                raise ValueError("label attribute %s not found" % scores)

        X2 = np.array(scores)

        XI = _an.XICOR(X1, X2)
        Z = np.sign(scale_matrix(X1).T @ scale_matrix(X2)) * XI["Z"]
        Enrichment = pd.DataFrame(Z, index=np.arange(1, X1.shape[1] + 1), columns=scores.columns)

        annotations = scores.columns
        idx = np.argmax(Z, axis=1)
        Label = annotations[idx]
        Label.index = np.arange(1, len(Label) + 1)
        Confidence = np.max(Z, axis=1)

    return Label, Confidence, Enrichment


def map_cell_scores(
    adata: AnnData,
    enrichment: pd.DataFrame,
    archetypes_key: Optional[str] = "H_unified",
    normalize: Optional[bool] = False,
):
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")

    enrichment_mat = np.array(enrichment)
    cell_scores_mat = adata.obsm[archetypes_key]

    if enrichment_mat.shape[0] != cell_scores_mat.shape[1]:
        raise ValueError(f"The number of rows in matrix `enrichment` ({enrichment.shape[0]}) must equal " f"the number of columns in matrix adata.obsm['{archetypes_key}'] ({cell_scores_mat.shape[1]}).")

    if normalize:
        enrichment_mat = double_normalize(enrichment_mat)

    cell_enrichment_mat = cell_scores_mat @ enrichment_mat
    Enrichment = pd.DataFrame(cell_enrichment_mat, index=adata.obs.index, columns=enrichment.columns)

    annotations = enrichment.columns
    idx = np.argmax(cell_enrichment_mat, axis=1)
    Label = annotations[idx]
    Label.index = np.arange(1, len(Label) + 1)
    Confidence = np.max(cell_enrichment_mat, axis=1)

    return Label, Confidence, Enrichment


def construct_backbone(
    adata: AnnData,
    density: Optional[float] = 0.5,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
):
    adata = adata.copy() if copy else adata

    footprint = adata.obsm["archetype_footprint"].T
    arch_features = sparse.csr_matrix(_an.normalize_mat(footprint, 1))
    arch_graph = build(data=arch_features, density=density)

    C = adata.obsm["C_unified"]
    coors2D = adata.obsm["ACTIONet2D"]
    coors3D = adata.obsm["ACTIONet2D"]
    denovo_colors = adata.obsm["denovo_color"]

    arch_coors2D = C.T @ coors2D
    arch_coors3D = C.T @ coors3D
    arch_denovo_colors = C.T @ denovo_colors

    backbone = {
        "graph": arch_graph,
        "coordinates": arch_coors2D,
        "coordinates_3D": arch_coors3D,
        "colors": arch_denovo_colors,
    }

    if return_raw:
        return backbone
    else:
        adata.uns["metadata"]["backbone"] = backbone

    return adata if copy else None
