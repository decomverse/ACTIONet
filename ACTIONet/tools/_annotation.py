from typing import Literal, Optional, Tuple

import numpy as np
from anndata import AnnData
from scipy import sparse

from . import _imputation as imputation
from . import _normalization as normalization

def annotate_archetypes_using_labels(
    adata: AnnData,
    label_key: Optional[str] = 'cell_type',
    archetypes_key: Optional[str] = 'H_unified',
) -> Tuple[list, np.ndarray, np.ndarray]:
    """
    """
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f'Did not find adata.obsm[\'{key}\'].')
    labels = adata.obs[label_key]
    profile = adata.obsm[archetypes_key].T
    if sparse.issparse(profile):
        profile = profile.toarray()

    # Compute enrichment using t-statistics
    unique_labels = labels.unique()
    Z = np.zeros((profile.shape[0], len(unique_labels)))
    for i, label in enumerate(unique_labels):
        mask = labels == label
        class_profile = profile[:, mask]
        null_profile = profile[:, ~mask]

        n_class = np.sum(mask)
        n_null = np.sum(~mask)

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
    significance_key: Optional[str] = 'unified_feature_specificity',
    n_iters: Optional[int] = 1000,
    seed: Optional[int] = 0,
) -> Tuple[list, np.ndarray, np.ndarray]:
    if significance_key not in adata.varm.keys():
        raise ValueError(
            f'Did not find adata.varm[\'{significance_key}\']. '
            'Please run pp.compute_archetype_feature_specificity() first.'
        )
    specificity = np.log1p(adata.varm[significance_key]).T
    np.nan_to_num(specificity, copy=False, nan=0.0)

    unique_genes = set([g for genes in marker_genes for g in genes])
    mask = adata.var.index.isin(unique_genes)
    specificity_genes = adata.var.index[mask]
    specificity_panel = specificity[:, mask]

    np.random.seed(seed)
    Z = np.empty((specificity_panel.shape[0], len(marker_genes)))
    names = names or [f'Celltype {i+1}' for i in range(len(marker_genes))]
    for i, (name, genes, ds) in enumerate(zip(names, marker_genes, directions)):
        mask = specificity_genes.isin(genes)
        A = specificity_panel[:, mask]
        sgn = np.array([
            -1 if direction == '-' else 1 for direction, gene in zip(ds, genes)
            if gene in specificity_genes
        ]).reshape(-1, 1)
        stat = (A @ sgn).flatten()

        rand_stats = np.empty((specificity_panel.shape[0], n_iters))
        for j in range(n_iters):
            rand_samples = np.random.choice(specificity_panel.shape[1], sum(gene in specificity_genes for gene in genes))
            rand_samples = np.array([2, 3])
            rand_A = specificity_panel[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(rand_stats, axis=1, ddof=1)
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    archetype_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    return archetype_labels, confidences, Z


def annotate_cells_using_markers(
    adata: AnnData,
    marker_genes: list,
    directions: list,
    names: Optional[list] = None,
    method: Optional[Literal['diffusion', 'archetype']] = 'diffusion',
    alpha: Optional[float] = 0.85,
    archetypes_key: Optional[str] = 'H_unified',
    significance: Optional[Literal['upper', 'lower']] = 'upper',
    n_threads: Optional[int] = 0,
    n_iters: Optional[int] = 1000,
    n_iters_diffusion: Optional[int] = 5,
) -> Tuple[list, np.ndarray, np.ndarray]:
    unique_genes = list(set([g for genes in marker_genes for g in genes]))

    if method == 'diffusion':
        adata_imputed = imputation.impute_genes_using_network(
            adata, unique_genes, alpha, n_threads, n_iters_diffusion
        )
    elif method == 'archetype':
        adata_imputed = imputation.impute_specific_genes_using_archetypes(
            adata, unique_genes, archetypes_key, significance
        )
    else:
        adata_imputed = adata

    Z = np.empty((adata_imputed.shape[0], len(marker_genes)))
    names = names or [f'Celltype {i+1}' for i in range(len(marker_genes))]
    for i, (name, genes, ds) in enumerate(zip(names, marker_genes, directions)):
        mask = adata_imputed.var.index.isin(genes)
        A = adata_imputed.X[:, mask]
        sgn = np.array([
            -1 if direction == '-' else 1 for direction, gene in zip(ds, genes)
            if gene in adata_imputed.var.index
        ]).reshape(-1, 1)
        stat = (A @ sgn).flatten()

        rand_stats = np.empty((adata_imputed.shape[0], n_iters))
        for j in range(n_iters):
            rand_samples = np.random.choice(adata_imputed.shape[1], sum(gene in adata_imputed.var.index for gene in genes))
            rand_A = adata_imputed.X[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(rand_stats, axis=1, ddof=1)
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    cell_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    return cell_labels, confidences, Z


def map_cell_scores_from_archetype_enrichment(
    adata: AnnData,
    enrichment: np.ndarray,
    archetypes_key: Optional[str] = 'H_unified',
    normalize: Optional[bool] = False
) -> np.ndarray:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f'Did not find adata.obsm[\'{archetypes_key}\'].')
    cell_scores_mat = adata.obsm[archetypes_key]

    if enrichment.shape[0] != cell_scores_mat.shape[1]:
        raise ValueError(
            f'The number of rows in matrix `enrichment` ({enrichment.shape[0]}) must equal '
            f'the number of columns in matrix adata.obsm[\'{archetypes_key}\'] ({cell_scores_mat.shape[1]}).'
        )

    if normalize:
        enrichment_scaled = normalization.double_normalize(enrichment)
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
    archetypes_key: Optional[str] = 'H_unified',
    significance_key: Optional[str] = 'unified_feature_specificity',
    n_iters: Optional[int] = 1000,
) -> Tuple[list, np.ndarray, np.ndarray]:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f'Did not find adata.obsm[\'{archetypes_key}\'].')

    names = names or [f'Celltype {i+1}' for i in range(len(marker_genes))]
    _, _, Z = annotate_archetypes_using_markers(
        adata, marker_genes, directions, names, significance_key, n_iters
    )
    cell_enrichment_mat = map_cell_scores_from_archetype_enrichment(
        adata, Z, archetypes_key, normalize=True
    )
    cell_labels = [names[i] for i in np.argmax(cell_enrichment_mat, axis=1)]
    confidences = np.max(Z, axis=1)
    return cell_labels, confidences, Z
