from typing import Optional, Union, Literal, Tuple

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from .. import tools as tl


def annotate_archetypes_using_labels(
        adata: AnnData,
        label_key: Optional[str] = "cell_types",
        archetypes_key: Optional[str] = "H_unified",
        ) -> Tuple[list, np.ndarray, np.ndarray]:
    """"""
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{key}'].")
    labels = adata.obs[label_key]
    profile = adata.obsm[archetypes_key].T

    # Compute enrichment using t-statistics
    unique_labels = labels.unique()
    Z = np.zeros((profile.shape[0], len(unique_labels)))
    for i, label in enumerate(unique_labels):
        mask = labels == label
        class_profile = profile[:, mask]
        null_profile = profile[:, ~mask]

        n_class = class_profile.shape[1]
        n_null = null_profile.shape[1]

        if n_class < 3 or n_null < 3:
            # Leave this column as zeros
            continue

        mu_class = np.mean(class_profile, axis=1)
        mu_null = np.mean(null_profile, axis=1)
        sigma_sq_class = np.var(class_profile, axis=1)
        sigma_sq_null = np.var(null_profile, axis=1)

        delta_mean = mu_class - mu_null
        t_statistic = delta_mean / np.sqrt(
                (sigma_sq_class / n_class) + (sigma_sq_null / n_null)
                )
        Z[:, i] = t_statistic

    archetype_labels = list(unique_labels[np.argmax(Z, axis=1)])
    confidences = np.max(Z, axis=1)
    return archetype_labels, confidences, Z


def annotate_archetypes_using_markers(
        adata: AnnData,
        marker_genes: list,
        directions: list,
        names: Optional[list] = None,
        archetypes_key: Optional[str] = "H_unified",
        significance: Optional[Literal["upper", "lower"]] = "upper",
        n_iters: Optional[int] = 1000,
        ) -> Tuple[list, np.ndarray, np.ndarray]:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")
    if f"{archetypes_key}_{significance}_significance" not in adata.varm.keys():
        raise ValueError(
                f"Did not find adata.varm['{archetypes_key}_{significance}_significance']. "
                "Please run pp.compute_archetype_feature_specificity() first."
                )
    specificity = np.log1p(
            adata.varm[f"{archetypes_key}_{significance}_significance"]
            ).T
    np.nan_to_num(specificity, copy=False, nan=0.0)

    unique_genes = set([g for genes in marker_genes for g in genes])
    mask = adata.var.index.isin(unique_genes)
    specificity_genes = adata.var.index[mask]
    specificity_panel = specificity[:, mask]

    Z = np.empty((specificity_panel.shape[0], len(marker_genes)))
    names = names or [f"Celltype {i + 1}" for i in range(len(marker_genes))]
    for i, (name, genes, ds) in enumerate(zip(names, marker_genes, directions)):
        mask = specificity_genes.isin(genes)
        A = specificity_panel[:, mask]
        sgn = np.array(
                [
                    -1 if direction == "-" else 1
                    for direction, gene in zip(ds, genes)
                    if gene in specificity_genes
                    ]
                ).reshape(-1, 1)
        stat = (A @ sgn).flatten()

        rand_stats = np.empty((specificity_panel.shape[0], n_iters))
        for j in range(n_iters):
            rand_samples = np.random.choice(specificity_panel.shape[1], np.sum(mask))
            rand_A = specificity_panel[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(
                rand_stats, axis=1, ddof=1
                )
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    archetype_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    return archetype_labels, confidences, Z


def map_cell_scores_from_archetype_enrichment(
        adata: AnnData,
        enrichment: np.ndarray,
        archetypes_key: Optional[str] = "H_unified",
        normalize: Optional[bool] = False,
        ) -> np.ndarray:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")
    cell_scores_mat = adata.obsm[archetypes_key]

    if enrichment.shape[0] != cell_scores_mat.shape[1]:
        raise ValueError(
                f"The number of rows in matrix `enrichment` ({enrichment.shape[0]}) must equal "
                f"the number of columns in matrix adata.obsm['{archetypes_key}'] ({cell_scores_mat.shape[1]})."
                )

    if normalize:
        enrichment_scaled = tl.double_normalize(enrichment)
    else:
        enrichment_scaled = enrichment.copy()
        enrichment_scaled[enrichment_scaled < 0] = 0
        if np.max(enrichment_scaled) > 50:
            enrichment_scaled = np.log1p(enrichment_scaled)

    cell_enrichment_mat = cell_scores_mat @ enrichment_scaled
    return cell_enrichment_mat


def annotate_cells_from_archetypes_using_markers(
        adata: AnnData,
        marker_genes: list,
        directions: list,
        names: Optional[list] = None,
        archetypes_key: Optional[str] = "H_unified",
        significance: Optional[Literal["upper", "lower"]] = "upper",
        n_iters: Optional[int] = 1000,
        ) -> Tuple[list, np.ndarray, np.ndarray]:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")

    names = names or [f"Celltype {i + 1}" for i in range(len(marker_genes))]
    _, _, Z = annotate_archetypes_using_markers(
            adata, marker_genes, directions, names, archetypes_key, significance, n_iters
            )
    cell_enrichment_mat = map_cell_scores_from_archetype_enrichment(
            adata, Z, archetypes_key, normalize=True
            )
    cell_labels = [names[i] for i in np.argmax(cell_enrichment_mat, axis=1)]
    confidences = np.max(Z, axis=1)
    return cell_labels, confidences, Z
