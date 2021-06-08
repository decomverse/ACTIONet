from typing import Optional
from typing_extensions import Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csc_matrix

import _ACTIONet as _an


def impute_genes_using_archetypes(
    adata: AnnData, genes: list, archetypes_key: Optional[str] = "H_unified"
) -> AnnData:

    """
    Impute expression of genes by interpolating over archetype profile

    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    genes
        List of genes to impute
    archetypes_key
        Key in `adata.obsm` that holds the archetype footprints

    Returns
    -------
        AnnData
            cells x genes
    """

    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")
    if f"{archetypes_key}_profile" not in adata.varm.keys():
        raise ValueError(
            f"Did not find adata.varm['{archetypes_key}_profile']. "
            "Please run pp.compute_archetype_feature_specificity() first."
        )

    genes = adata.obs.index.intersection(genes)
    Z = adata[:, genes].varm[f"{archetypes_key}_profile"]
    H = adata.obsm[archetypes_key].T
    return AnnData(
        X=(Z @ H).T,
        obs=pd.DataFrame(index=adata.obs.index),
        var=pd.DataFrame(index=genes),
    )


def impute_specific_genes_using_archetypes(
    adata: AnnData,
    genes: list,
    archetypes_key: Optional[str] = "H_unified",
    significance: Optional[Literal["upper", "lower"]] = "upper",
) -> AnnData:

    """
    Impute expression of genes by interpolating over archetype profile

    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    genes
        List of genes to impute
    archetypes_key:
        Key in `adata.obsm` that holds the archetype footprints
    significance:
        Whether to use upper or lower significant genes.

    Returns
    -------
        AnnData
            cells x genes
    """
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")
    if f"{archetypes_key}_{significance}_significance" not in adata.varm.keys():
        raise ValueError(
            f"Did not find adata.varm['{archetypes_key}_{significance}_significance']. "
            "Please run pp.compute_archetype_feature_specificity() first."
        )
    genes = adata.obs.index.intersection(genes)
    Z = np.log1p(adata[:, genes].varm[f"{archetypes_key}_{significance}_significance"])
    H = adata.obsm[archetypes_key].T
    return AnnData(
        X=(Z @ H).T,
        obs=pd.DataFrame(index=adata.obs.index),
        var=pd.DataFrame(index=genes),
    )


def impute_genes_using_network(
    adata: AnnData,
    genes: list,
    alpha: Optional[float] = 0.85,
    thread_no: Optional[int] = 0,
    n_iters: Optional[int] = 5,
) -> AnnData:
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            f"Did not find adata.obsp['ACTIONet']. "
            "Please run nt.build_network() first."
        )

    genes = adata.var.index.intersection(genes)
    mask = adata.var.index.isin(genes)
    if (np.sum(mask)) > 0:
        U = adata.X[:, mask].copy()
        U[U < 0] = 0
        cs = np.sum(U, axis=0)
        U = U / cs
        U = U[:, cs > 0]
        gg = genes[cs > 0]
    else:
        U = adata.X[:, mask].copy()
        U = U / np.sum(U)
        gg = genes

    # Network diffusion
    G = adata.obsp["ACTIONet"]
    imputed = _an.compute_network_diffusion(G, csc_matrix(U), thread_no, alpha, n_iters)
    np.nan_to_num(imputed, copy=False, nan=0.0)

    # Rescale the baseline expression of each gene
    rescaled = np.zeros(imputed.shape)
    for i in range(imputed.shape[1]):
        x = U[:, i]
        y = imputed[:, i]

        # In the R implementation, quantile(x, 1) is used, which just
        # finds the maximum.
        x_Q = np.max(x)
        y_Q = np.max(y)

        if y_Q == 0:
            # Leave this column as zeros
            continue

        y = y * x_Q / y_Q
        y[y > x_Q] = x_Q
        rescaled[:, i] = y

    return AnnData(
        X=rescaled,
        obs=pd.DataFrame(index=adata.obs.index),
        var=pd.DataFrame(index=genes),
    )
