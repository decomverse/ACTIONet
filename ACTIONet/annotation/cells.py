from typing import Optional, Tuple, Literal

import numpy as np
from anndata import AnnData

from ..tools import broken_shit


def _annotate_cells_using_markers(
        adata: AnnData,
        marker_genes: list,
        directions: list,
        names: Optional[list] = None,
        method: Optional[Literal["diffusion", "archetype"]] = "diffusion",
        alpha: Optional[float] = 0.85,
        archetypes_key: Optional[str] = "H_unified",
        significance: Optional[Literal["upper", "lower"]] = "upper",
        thread_no: Optional[int] = 0,
        n_iters: Optional[int] = 1000,
        n_iters_diffusion: Optional[int] = 5
        ) -> Tuple[list, np.ndarray, np.ndarray]:

    unique_genes = list(set([g for genes in marker_genes for g in genes]))

    if method == "diffusion":
        adata_imputed = imputation.impute_genes_using_network(
                adata, unique_genes, alpha, thread_no, n_iters_diffusion
                )
    elif method == "archetype":
        adata_imputed = imputation.impute_specific_genes_using_archetypes(
                adata, unique_genes, archetypes_key, significance
                )
    else:
        adata_imputed = adata

    Z = np.empty((adata_imputed.shape[0], len(marker_genes)))
    names = names or [f"Celltype {i + 1}" for i in range(len(marker_genes))]
    for i, (name, genes, ds) in enumerate(zip(names, marker_genes, directions)):
        mask = adata_imputed.var.index.isin(genes)
        A = adata_imputed.X[:, mask]
        sgn = np.array(
                [
                    -1 if direction == "-" else 1
                    for direction, gene in zip(ds, genes)
                    if gene in adata_imputed.var.index
                    ]
                ).reshape(-1, 1)
        stat = (A @ sgn).flatten()

        rand_stats = np.empty((adata_imputed.shape[0], n_iters))
        for j in range(n_iters):
            rand_samples = np.random.choice(adata_imputed.shape[1], np.sum(mask))
            rand_A = adata_imputed.X[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(
                rand_stats, axis=1, ddof=1
                )
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    cell_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    return cell_labels, confidences, Z
