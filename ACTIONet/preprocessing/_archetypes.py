from typing import Optional, Union
from typing_extensions import Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy import sparse

import _ACTIONet as _an


def prune_archetypes(
    adata: AnnData,
    C_trace: list,
    H_trace: list,
    min_specificity_z_threshold: Optional[float] = -3,
    min_cells: Optional[int] = 2,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """\
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    adata:
        Current AnnData object storing the ACTIONet results
    C_trace, H_trace:
        Output of pp.ACTION()
    min_specificity_z_threshold:
        Controls level of prunning for non-specific archetypes
        (larger values remove more archetypes)
    min_cells:
        Minimum number of influential cells for each archetype
        to be considered nontrivial
    copy
        Determines whether a copy of `adata` is returned.
    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['ACTION_C_stacked']`
        `.obsm['ACTION_H_stacked']`
        `.uns['ACTION']['archetypes']['pruned']`
    """
    if "ACTION" not in adata.uns.keys():
        raise ValueError(
            "Did not find adata.uns['ACTION']. " "Please run pp.ACTION() first."
        )

    pruned = _an.prune_archetypes(
        C_trace, H_trace, min_specificity_z_threshold, min_cells
    )

    adata = adata.copy() if copy else adata
    adata.obsm["ACTION_C_stacked"] = sparse.csc_matrix(pruned["C_stacked"])
    adata.obsm["ACTION_H_stacked"] = sparse.csc_matrix(pruned["H_stacked"].T)
    adata.uns["ACTION"].setdefault("archetypes", {}).update(
        {"pruned": {"selected_archetypes": pruned["selected_archs"]}}
    )

    return adata if copy else None


def unify_archetypes(
    adata: AnnData,
    alpha: Optional[float] =  0.99, 
    outlier_threshold: Optional[float] = 2, 
    sim_threshold: Optional[float] = 0, 
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> AnnData:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    adata:
        Current AnnData object storing the ACTIONet results

    alpha: 
        Diffusion parameter to impute archetype foorprints ([0-1), default: 0.99)
                                                             
    outlier_threshold: 
        Coreness threshold to filter noisy archetypes(<= 0, default: 2)
        
    sim_threshold: 
        Similarity threshold to group similar archetypes (default: 0)

    n_threads:
        Number of threads (default: 0 [all])
        
    copy
        Determines whether a copy of `adata` is returned.
    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['ACTION_C_unified']`
        `.obsm['ACTION_H_unified']`
        `.uns['ACTION']['archetypes']['unified']`
    """
    # Check for ACTION_S_r and ACTION_H_stacked
    if "ACTION_S_r" not in adata.obsm.keys():
        raise ValueError(
            "Did not find adata.obsm['ACTION_S_r']. "
            "Please run pp.reduce_kernel() first."
        )
    if (
        "ACTION_C_stacked" not in adata.obsm.keys()
        or "ACTION_H_stacked" not in adata.obsm.keys()
    ):
        raise ValueError(
            "Did not find adata.obsm['ACTION_C_stacked'] or adata.obsm['ACTION_H_stacked']. "
            "Please run pp.prune_archetypes() first."
        )

    adata = adata.copy() if copy else adata
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            "Did not find adata.obsp['ACTIONet']. "
            "Please run nt.build_network() first."
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp["ACTIONet"]
    S_r = adata.obsm["ACTION_S_r"].T
    C = adata.obsm["ACTION_C_stacked"]
    if sparse.issparse(C):
        C = C.toarray()
    H = adata.obsm["ACTION_H_stacked"].T
    if sparse.issparse(H):
        H = H.toarray()
    unified = _an.unify_archetypes(
        G,
        S_r,
        C,
        alpha,
        outlier_threshold,
        sim_threshold,
        n_threads
    )

    adata.obsm["ACTION_C_unified"] = sparse.csc_matrix(unified["C_unified"])
    adata.obsm["ACTION_H_unified"] = sparse.csc_matrix(unified["H_unified"].T)
    adata.uns["ACTION"].setdefault("archetypes", {}).update(
        {"unified": {"selected_archetypes": unified["selected_archetypes"]}}
    )

    groups = unified["assigned_archetypes"]
    adata.obs["ACTION"] = pd.Categorical(
        values=groups.astype("U"),
        categories=natsorted(map(str, np.unique(groups))),
    )

    return adata if copy else None
