from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy import sparse

import _ACTIONet as _an


def prune_archetypes(
        C_trace: list,
        H_trace: list,
        adata: Optional[AnnData] = None,
        min_specificity_z_threshold: Optional[float] = -3,
        min_cells: Optional[int] = 2,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False
) -> Union[dict, AnnData, None]:
    """\
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    C_trace, H_trace:
        Output of pp.ACTION()
    adata
        Optional AnnData object in which to store output of 'prune_archetypes()'
    min_specificity_z_threshold:
        Controls level of prunning for non-specific archetypes
        (larger values remove more archetypes)
    min_cells:
        Minimum number of influential cells for each archetype
        to be considered nontrivial
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'return_raw=True' and 'data' is AnnData, return raw output of 'prune_archetypes()'.
    Returns
    -------

    pruned : dict
        dict containing 'C_stacked' and 'H_stacked' matrices

    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm['C_stacked']`
        `.obsm['H_stacked']`
        `.uns['obsm_annot']['C_stacked']['type']`
        `.uns['obsm_annot']['H_stacked']['type']`
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
        else:
            raise ValueError("'adata' is not an AnnData object.")

    pruned = _an.prune_archetypes(C_trace, H_trace, min_specificity_z_threshold, min_cells)

    pruned["C_stacked"] = sparse.csc_matrix(pruned["C_stacked"])
    pruned["H_stacked"] = sparse.csc_matrix(pruned["H_stacked"].T)

    if return_raw or adata is None:
        return pruned
    else:
        adata.obsm["C_stacked"] = pruned["C_stacked"]
        adata.obsm["H_stacked"] = pruned["H_stacked"]
        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "C_stacked": {"type": np.array([b'internal'], dtype=object)},
                "H_stacked": {"type": np.array([b'internal'], dtype=object)}
            })
        return adata if copy else None


def unify_archetypes(
        adata: Optional[AnnData] = None,
        G: Union[np.ndarray, sparse.spmatrix] = None,
        S_r: Union[np.ndarray, sparse.spmatrix] = None,
        C_stacked: Union[np.ndarray, sparse.spmatrix] = None,
        reduction_name: Optional[str] = "ACTION",
        net_name: Optional[str] = "ACTIONet",
        alpha_val: Optional[float] = 0.99,
        outlier_threshold: Optional[float] = 2,
        sim_threshold: Optional[float] = 0,
        thread_no: Optional[int] = 0,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False,
) -> AnnData:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    adata:
        AnnData object containing 'reduction_name' in '.obsm' and 'net_name' in '.obsp'.
    G:
        Adjacency matrix to use for unification.
        Required if 'adata=None'.
    S_r:
        Reduced representation matrix to use for unification.
        Required if 'adata=None'.
    C_stacked:
        Matrix containing output 'C_stacked' of 'prune_archetypes()' to use for unification.
        Required if 'adata=None', otherwise retrieved from 'adata.obsm["C_stacked"]'
    reduction_name:
        Key of 'adata.obms' containing reduced matrix to use for 'S_r' in 'unify_archetypes()' (default="ACTION").
        Ignored if 'adata=None'.
    net_name:
        Key of 'adata.obmp' containing adjacency matrix to use for 'G' in 'unify_archetypes()' (default="ACTIONet").
        Ignored if 'adata=None'.
    alpha_val:
        Diffusion parameter to impute archetype foorprints ([0-1), default: 0.99)
    outlier_threshold:
        Coreness threshold to filter noisy archetypes(<= 0, default: 2)
    sim_threshold: 
        Similarity threshold to group similar archetypes (default: 0)
    thread_no:
        Number of threads. Defaults to number of threads available - 2.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return dict of raw 'unify_archetypes()' output instead of storing to 'adata'.
    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['C_unified']`
        `.obsm['H_unified']`
        `.uns['ACTION']['archetypes']['unified']`

    unified : dict
        If 'adata=None' or 'return_raw=True', a dict containing 'C_unified' and 'H_unified' matrices
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            S_r = S_r if S_r is not None else adata.obsm[reduction_name]
            G = G if G is not None else adata.obsp[net_name]
            C_stacked = C_stacked if C_stacked is not None else adata.obsm['C_stacked']
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if G is None or S_r is None or C_stacked is None:
            raise ValueError("'G' and 'S_r' cannot be NoneType if 'adata=None'.")
        if not isinstance(G, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(S_r, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'S_r' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(C_stacked, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'C_stacked' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    S_r = S_r.T.astype(dtype=np.float64)
    C_stacked = C_stacked.astype(dtype=np.float64)

    if sparse.issparse(C_stacked):
        C_stacked = C_stacked.toarray()

    unified = _an.unify_archetypes(G, S_r, C_stacked, alpha_val, outlier_threshold, sim_threshold, thread_no)

    unified["C_unified"] = sparse.csc_matrix(unified["C_unified"])
    unified["H_unified"] = sparse.csc_matrix(unified["H_unified"].T)

    if return_raw or adata is None:
        return unified
    else:
        adata.obsm["C_unified"] = unified["C_unified"]
        adata.obsm["H_unified"] = unified["H_unified"]

        # adata.uns.setdefault("metadata", {}).update(
        #     {"unified": {"selected_archetypes": unified["selected_archetypes"]}}
        # )
        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "C_unified": {"type": np.array([b'internal'], dtype=object)},
                "H_unified": {"type": np.array([b'internal'], dtype=object)},
            })

        groups = unified["assigned_archetype"]
        adata.obs["assigned_archetype"] = pd.Categorical(
            values=groups.astype(int),
            categories=natsorted(map(int, np.unique(groups))),
        )

        return adata if copy else None

