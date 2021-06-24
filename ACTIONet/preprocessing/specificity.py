from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def _compute_archetype_specificity(S, H):
    H = np.array(H, dtype=np.float64)
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        return _an.compute_archetype_feature_specificity(S, H)
    return _an.compute_archetype_feature_specificity_full(S, H)


def _compute_cluster_specificity(S, assignments):
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        return _an.compute_cluster_feature_specificity(S, assignments)
    return _an.compute_cluster_feature_specificity_full(S, assignments)


def compute_archetype_feature_specificity(
    adata: Optional[AnnData] = None,
    S: Union[np.ndarray, sparse.spmatrix] = None,
    H: Union[np.ndarray, sparse.spmatrix] = None,
    layer_name: Optional[str] = None,
    footprint_key: Optional[str] = "archetype_footprint",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False
) -> Union[AnnData, dict, None]:

    """Computes Feature (i.e., gene) specificity of archetypes \
    Uses Archetype footprints to estimate markers (soft clustering) 
    :param adata: AnnData object possibly containing '.layers["layer_name"]' and '.obsm["footprint_key"]'.
    :param S: `n_obs` × `n_vars` gene expression matrix. \
        Required if 'adata=None', otherwise retrieved from '.layers["layer_name"]' or '.X' if 'layer_name=None'.
    :param H: Matrix-like object containing archetype footprint \
        Required if 'adata=None', otherwise retrieved from '.obsm["footprint_key"]'
    :param layer_name: Key of 'layers' to retrieve gene expression matrix.
    :param footprint_key:Key in `adata.obsm` that holds the archetype footprint.
    :param copy: If 'adata' is given, return a copy instead of writing to `adata`
    :param return_raw: If 'adata' is given, return dict of 'compute_archetype_feature_specificity()' output instead of storing to 'adata'.
    ...
    :return adata : anndata.AnnData \
        if `copy=True` returns None or else adds fields to `adata` \
        `.varm["unified_feature_profile"]` \
        `.varm["unified_feature_specificity"]`
    :return specificity: dict \
        If 'adata=None' or 'return_raw=True', returns dict containing 'unified_feature_profile' and 'unified_feature_specificity' matrices.
    """
    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            if layer_name is not None:
                S = S.T if S is not None else adata.layers[layer_name].T
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

    specificity = _compute_archetype_specificity(S, H)

    if return_raw or adata is None:
        return specificity
    else:
        adata.varm["unified_feature_profile"] = specificity["archetypes"]
        adata.varm["unified_feature_specificity"] = specificity["upper_significance"]
        # adata.varm["lower_significance"] = specificity["lower_significance"]

        adata.uns.setdefault("varm_annot", {}).update(
            {
                "unified_feature_profile": {"type": np.array([b'internal'], dtype=object)},
                "unified_feature_specificity": {"type": np.array([b'reduction'], dtype=object)},
            })

        return adata if copy else None


def compute_cluster_feature_specificity(
    adata: AnnData, cluster_key: Optional[str] = "leiden", copy: Optional[bool] = False
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of clusters

    Uses cluster membership vector to estimate markers (disjoint clustering)
    Parameters
    ----------
    adata
        Current AnnData object storing the ACTIONet results
    cluster_key
        Key in `adata.obs` that holds the clustering variable
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.varm[f'{cluster_key}_profile']`
        `.varm[f'{cluster_key}_upper_significance']`
        `.varm[f'{cluster_key}_lower_significance']`
    """
    if cluster_key not in adata.obs.keys():
        raise ValueError(f"Did not find adata.obs['{cluster_key}'].")
    adata = adata.copy() if copy else adata

    S = adata.X.T
    assignments = adata.obs[cluster_key]
    if isinstance(assignments, pd.Series):
        assignments = pd.factorize(assignments)[0]

    specificity = _compute_cluster_specificity(S, assignments)

    adata.varm[f"{cluster_key}_profile"] = specificity["archetypes"]
    adata.varm[f"{cluster_key}_upper_significance"] = specificity["upper_significance"]
    adata.varm[f"{cluster_key}_lower_significance"] = specificity["lower_significance"]

    return adata if copy else None
