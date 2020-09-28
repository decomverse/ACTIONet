from typing import Literal, Optional

import numpy as np
from anndata import AnnData

from . import _imputation as imputation

def annotate_archetypes_using_labels(
    adata: AnnData,
    label_key: Optional[str] = 'cell_types',
    archetypes_key: Optional[str] = 'ACTION_H_unified',
    copy: Optional[bool] = False,
) -> AnnData:
    """
    """
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f'Did not find adata.obsm[\'{key}\'].')
    adata = adata.copy() if copy else adata
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

    archetype_labels = unique_labels[np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    adata.uns['ACTION']['archetypes']['label_annotation'] = {
        'archetype_labels': archetype_labels,
        'confidences': confidences,
        'enrichment': Z,
        'label_key': label_key,
        'archetypes_key': archetypes_key,
    }

    return adata if copy else None


def annotate_archetypes_using_markers(
    adata: AnnData,
    marker_genes: list,
    directions: list,
    names: Optional[list] = None,
    archetypes_key: Optional[str] = 'ACTION_H_unified',
    significance: Optional[Literal['upper', 'lower']] = 'upper',
    n_iters: Optional[int] = 1000,
    copy: Optional[bool] = False,
) -> AnnData:
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f'Did not find adata.obsm[\'{archetypes_key}\'].')
    if f'{archetypes_key}_{significance}_significance' not in adata.varm.keys():
        raise ValueError(
            f'Did not find adata.varm[\'{archetypes_key}_{significance}_significance\']. '
            'Please run pp.compute_archetype_feature_specificity() first.'
        )
    adata = adata.copy() if copy else adata
    specificity = np.log1p(adata.varm[f'{archetypes_key}_{significance}_significance']).T
    np.nan_to_num(specificity, copy=False, nan=0.0)

    unique_genes = set([g for genes in marker_genes for g in genes])
    mask = adata.var.index.isin(unique_genes)
    specificity_genes = adata.var.index[mask]
    specificity_panel = specificity[:, mask]

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
            rand_samples = np.random.choice(specificity_panel.shape[1], np.sum(mask))
            rand_A = specificity_panel[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(rand_stats, axis=1, ddof=1)
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    archetype_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    adata.uns['ACTION']['archetypes']['marker_annotation'] = {
        'archetype_labels': archetype_labels,
        'confidences': confidences,
        'enrichment': Z,
        'params': {
            'marker_genes': marker_genes,
            'directions': directions,
            'names': names,
            'archetypes_key': archetypes_key,
        }
    }

    return adata if copy else None


def annotate_cells_using_markers(
    adata: AnnData,
    marker_genes: list,
    directions: list,
    names: Optional[list] = None,
    method: Optional[Literal['diffusion', 'archetype']] = 'diffusion',
    alpha: Optional[float] = 0.85,
    archetypes_key: Optional[str] = 'ACTION_H_unified',
    significance: Optional[Literal['upper', 'lower']] = 'upper',
    n_threads: Optional[int] = 0,
    n_iters: Optional[int] = 1000,
    n_iters_diffusion: Optional[int] = 5,
    copy: Optional[bool] = False,
) -> AnnData:
    unique_genes = list(set([g for genes in marker_genes for g in genes]))
    adata = adata.copy() if copy else adata

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
            rand_samples = np.random.choice(adata_imputed.shape[1], np.sum(mask))
            rand_A = adata_imputed.X[:, rand_samples]
            rand_stat = (rand_A @ sgn).flatten()
            rand_stats[:, j] = rand_stat

        cell_zscores = (stat - np.mean(rand_stats, axis=1)) / np.std(rand_stats, axis=1, ddof=1)
        Z[:, i] = cell_zscores
    np.nan_to_num(Z, copy=False, nan=0.0)

    cell_labels = [names[i] for i in np.argmax(Z, axis=1)]
    confidences = np.max(Z, axis=1)
    adata.uns['ACTION'].setdefault('cells', {})['marker_annotation'] = {
        'cell_labels': cell_labels,
        'confidences': confidences,
        'enrichment': Z,
        'params': {
            'marker_genes': marker_genes,
            'directions': directions,
            'names': names,
            'method': method,
            'archetypes_key': archetypes_key,
        }
    }
    return adata if copy else None

def map_cell_scores_from_archetype_enrichment():
    pass

def annotate_cells_from_archetypes_using_markers():
    pass
