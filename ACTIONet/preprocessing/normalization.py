import warnings
from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def normalize_matrix(
    X: Union[np.ndarray, sparse.spmatrix],
    top_features_frac: float = 1.0,
    scale_factor: Union[str, float, int, np.ndarray, None] = "median",
    transformation: Union[str, None] = "log",
    anchor_features: Union[np.ndarray, None] = None,
) -> Union[np.ndarray, sparse.spmatrix]:

    X = X.astype(dtype=np.float64)

    # Which features (i.e. genes) should we use to compute library sizes?
    if anchor_features is not None:
        lib_sizes = np.array(np.mean(X[:, anchor_features], axis=1))
    else:
        if top_features_frac < 1.0:
            universality = np.array(np.mean(X > 0, axis=0))
            selected_features = np.flatnonzero(universality > (1 - top_features_frac))
            lib_sizes = np.array(np.mean(X[:, selected_features], axis=1))
        else:
            lib_sizes = np.array(np.mean(X, axis=1))

    # Note: mean as opposed to sum

    # Normalize library sizes
    if isinstance(X, sparse.spmatrix):
        X_scaled = X.multiply(1 / lib_sizes)
    else:
        X_scaled == X / lib_sizes

    # scale normalized columns
    if scale_factor == "median":
        kappa = np.median(np.array(np.sum(X, axis=1) / np.sum(X_scaled, axis=1)))
        X_scaled_norm = X_scaled * kappa
    elif isinstance(scale_factor, (int, float)):
        X_scaled_norm = X_scaled * scale_factor
    elif isinstance(scale_factor, np.ndarray):
        if sparse.issparse(X_scaled):
            X_scaled_norm = X_scaled.multiply(scale_factor)
        else:
            X_scaled_norm = X_scaled / scale_factor

    # For compatibility with C
    if sparse.issparse(X_scaled_norm):
        X_scaled_norm = sparse.csc_matrix(X_scaled_norm)

    # Post-transformation
    if transformation == "log":
        X_scaled_norm_trans = np.log1p(X_scaled_norm)
    elif transformation == "tukey":
        if sparse.issparse(X_scaled_norm):
            nnz_idx = X_scaled_norm.nonzero()
            ii = nnz_idx[0]
            jj = nnz_idx[1]
            vv = X_scaled_norm[ii, jj]
            vv_transformed = np.sqrt(vv) + np.sqrt(1 + vv)
            X_scaled_norm[ii, jj] = vv_transformed
        else:
            X_scaled_norm[X_scaled_norm < 0] = 0
            vv = X_scaled_norm[X_scaled_norm != 0]
            vv_transformed = np.sqrt(vv) + np.sqrt(1 + vv)
            X_scaled_norm[X_scaled_norm != 0] = vv_transformed
    elif transformation == "lsi":
        if sparse.issparse(X_scaled_norm):
            X_scaled_norm_trans = _an.LSI(X_scaled_norm)
        else:
            X_scaled_norm_sp = sparse.csc_matrix(X_scaled_norm)
            X_scaled_norm_trans = _an.LSI(X_scaled_norm_sp).toarray()
    else:
        X_scaled_norm_trans = X_scaled_norm

    return X_scaled_norm_trans


def normalize(
    adata: AnnData,
    layer_key: Optional[str] = None,
    layer_key_out: Optional[str] = None,
    top_features_frac: float = 1.0,
    scale_factor: Union[str, float, int, np.ndarray, None] = "median",
    transformation: Union[str, None] = "log",
    anchor_features: Union[np.ndarray, None] = None,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata

    if "norm_method" in adata.uns["metadata"].keys():  # Already normalized? leave it alone!
        # return adata if copy else None
        warnings.warn("AnnData object is prenormalized. Please make sure to use the right assay.")

    if layer_key is None and "input_assay" in adata.uns["metadata"].keys():
        layer_key = adata.uns["metadata"]["input_assay"]

    if layer_key is not None:
        if layer_key not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_key + "']. ")
        S = adata.layers[layer_key]
    else:
        S = adata.X

    if sparse.issparse(S):
        UE = set(S.data)
    else:
        UE = set(S.flatten())

    nonint_count = len(UE.difference(set(np.arange(0, max(UE) + 1))))
    if 0 < nonint_count:
        warnings.warn("Input [count] assay has non-integer values, which looks like a normalized matrix. Please make sure to use the right assay.")

    S = normalize_matrix(
        S,
        anchor_features=anchor_features,
        top_features_frac=top_features_frac,
        scale_factor=scale_factor,
        transformation=transformation,
    )

    adata.uns["metadata"]["norm_method"] = "default_top%.2f_%s" % (
        top_features_frac,
        transformation,
    )

    if layer_key_out is not None:
        adata.uns["metadata"]["default_assay"] = layer_key_out
        adata.layers[layer_key_out] = S
    else:
        adata.uns["metadata"]["default_assay"] = None
        adata.X = S

    return adata if copy else None
