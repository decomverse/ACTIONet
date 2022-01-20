from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from .. import tools as tl


def __compute_archetype_specificity(S, H, thread_no=0):
    H = np.array(H, dtype=np.float64)
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        out = _an.compute_archetype_feature_specificity(S, H, thread_no)
    else:
        S = np.array(S, dtype=np.float64)
        out = _an.compute_archetype_feature_specificity_full(S, H, thread_no)

    return out


def __compute_cluster_specificity(S, sample_assignments, thread_no=0):
    if sparse.issparse(S):
        S = S.tocsc().astype(dtype=np.float64)
        out = _an.compute_cluster_feature_specificity(S, sample_assignments, thread_no)
    else:
        S = np.array(S, dtype=np.float64)
        out = _an.compute_cluster_feature_specificity_full(S, sample_assignments, thread_no)

    return out


def compute_archetype_feature_specificity(
        adata: Optional[AnnData] = None,
        S: Union[np.ndarray, sparse.spmatrix] = None,
        H: Union[np.ndarray, sparse.spmatrix] = None,
        layer_key: Optional[str] = None,
        footprint_key: Optional[str] = "archetype_footprint",
        thread_no: Optional[int] = 0,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False
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

    specificity = __compute_archetype_specificity(S=S, H=H, thread_no=thread_no)

    if return_raw or adata is None:
        return specificity
    else:
        adata.varm["unified_feature_profile"] = specificity["archetypes"]
        adata.varm["unified_feature_specificity"] = specificity["upper_significance"]

        adata.uns.setdefault("varm_annot", {}).update(
                {
                    "unified_feature_profile": {"type": np.array([b'internal'], dtype=object)},
                    "unified_feature_specificity": {"type": np.array([b'reduction'], dtype=object)},
                    }
                )

        return adata if copy else None


def compute_cluster_feature_specificity(
        adata: AnnData,
        cluster_key: Optional[str] = None,
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
    adata
        Current AnnData object storing the ACTIONet results.
    cluster_key
        Key in `adata.obs` that holds the clustering variable.
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

    if cluster_key not in adata.obs.keys():
        raise ValueError(f"'{cluster_key}' not in adata.obs_keys()")
    adata = adata.copy() if copy else adata

    if layer_key is not None:
        S = adata.layers[layer_key]
    else:
        S = adata.X

    S = S.T.astype(dtype=np.float64)

    idx_clust = tl.get_data_or_split(adata=adata, attr=cluster_key, to_return="levels")

    specificity_out = __compute_cluster_specificity(S=S, sample_assignments=idx_clust["index"], thread_no=thread_no)

    if return_raw:
        return specificity_out
    else:
        adata.varm[f"{output_prefix}_feature_specificity"] = specificity_out["upper_significance"]
        adata.uns.setdefault("varm_annot", {}).update(
                {
                    f"{output_prefix}_feature_specificity": {"type": np.array([b'reduction'], dtype=object)},
                    }
                )

        return adata if copy else None
