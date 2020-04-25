import _ACTIONet as an
import numpy as np
import pandas as pd
from scipy.sparse import issparse, spmatrix
from typing import Optional, Union
from natsort import natsorted
from anndata import AnnData



def reduce_kernel(
    adata: Union[AnnData, np.ndarray, spmatrix],
    n_comps: Optional[int] = 50,
    svd_solver: Optional[str] = 'Halko',
    random_state: Optional[int] = 0,
    return_info: bool = False,
    use_highly_variable: Optional[bool] = None,
    dtype: str = 'float32',
    copy: bool = False
) -> [AnnData, np.ndarray, spmatrix]:
    """\
    Kernel Reduction Method [Mohammadi2020].

    Computes SVD-reduced form of the kernel matrix.

    Parameters
    ----------
    adata
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute. Defaults to 50, or 1 - minimum
        dimension size of selected representation.
    svd_solver
        SVD solver to use:
        `'Randomized'` (the default)
          Randomized SVD algorithm from: "Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions"
    seed
        Random seed        
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.        
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.        
    dtype
        Numpy data type string to which to convert the result.        
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    ACTION_S_r : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray`
        If `data` is array-like and `return_info=False` was passed,
        this function only returns `ACTION_S_r`…
    adata : anndata.AnnData
        …otherwise if `copy=True` it returns or else adds fields to `adata`:

        `.obsm['ACTION_S_r']`
             Scaled right singular vectors (reduced cell representations)
        `.varm['ACTION_V']`
             Left singular vectors (signifying gene modules)
        `.uns['ACTION_reduction']['lambda']`
             sigma_sq / (n-1)
        `.uns['ACTION_reduction']['explained_var']`
             Explained variance.
    """
    data_is_AnnData = isinstance(adata, AnnData)
    if data_is_AnnData:
        adata = adata.copy() if copy else adata
    else:
        adata = AnnData(adata)

    if use_highly_variable is True and 'highly_variable' not in adata.var.keys():
        raise ValueError(
            'Did not find adata.var[\'highly_variable\']. '
            'Either your data already only consists of highly-variable genes '
            'or consider running `pp.highly_variable_genes` first.'
        )
    if use_highly_variable is None:
        use_highly_variable = True if 'highly_variable' in adata.var.keys() else False
    adata_comp = (
        adata[:, adata.var['highly_variable']] if use_highly_variable else adata
    )

    X = adata_comp.X.T

    if (issparse(X)):    
        print("Running reduction in sparse mode ...")    
        reduction_out = an.reduce_kernel_sparse(X)
    else:
        print("Running reduction in dense mode ...")    
        reduction_out = an.reduce_kernel(X)
        
    ACTION_S_r = reduction_out['S_r'].T
    if ACTION_S_r.dtype.descr != np.dtype(dtype).descr:
        ACTION_S_r = ACTION_S_r.astype(dtype)

    print("Done.")    

    if data_is_AnnData:
        adata.obsm['ACTION_S_r'] = ACTION_S_r
        adata.uns['ACTIONet_reduction'] = {}
        adata.uns['ACTIONet_reduction']['params'] = {
            'genes_subset': adata.var['highly_variable'] if use_highly_variable else None
        }
        if use_highly_variable:
            adata.varm['ACTION_V'] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm['ACTION_V'][adata.var['highly_variable']] = reduction_out['V']
        else:
            adata.varm['ACTION_V'] = reduction_out['V']
            
        adata.uns['ACTIONet_reduction']['lambda'] = reduction_out['lambda']
        adata.uns['ACTIONet_reduction']['explained_var'] = reduction_out['explained_var']

        return adata if copy else None
    else:
        if return_info:
            return (
                ACTION_S_r,
                reduction_out['V'].T,
                reduction_out['lambda'],
                reduction_out['explained_var']
            )
        else:
            return ACTION_S_r


def run_simplex_regression(
    A: np.ndarray,
    B: np.ndarray
) -> np.ndarray:
    """\
    Simplex-Constrained Regression (AX-B).

    Solves the linear regression problem, subject to coefficients being positive and sum to one.

    Parameters
    ----------
    A:
        Matrix of independent variables (design matrix)
    B:
        Matrix of dependent variables (response variable)

    Returns
    -------
    X:
        Coefficient matrix
    """
    
    X = an.run_simplex_regresion(A, B)        
    return X


def run_SPA(
    A: np.ndarray,
    k: int = 50
) -> np.ndarray:
    """\
    Successive Projection Algorithm (SPA).

    Runs SPA algorithm to solve separable NMF problem.

    Parameters
    ----------
    A:
        Matrix matrix
    k:
        Number of columns to select

    Returns
    -------
    selected_columns:
        Index of k selected columns
    norms:
        Residual norm of selected columns    
    """
    
    SPA_out = an.run_SPA(A, k)            
    return(SPA_out["selected_columns"], SPA_out["norms"])



def run_AA(
    A: np.ndarray,
    k: int,
    max_iter: Optional[int] = 50,
    min_delta: Optional[float] = 0.01
) -> np.ndarray:
    """\
    Archetypal Analysis (AA)

    Runs SPA algorithm to solve separable NMF problem.

    Parameters
    ----------
    A:
        Input matrix
    W0:
        Initial estimate of archetypes with k (# archetype) columns
    max_it, min_delta:
        Define stopping conditions

    Returns
    -------
    C:
        Convex matrix of archetype coefficients (#observations x # archetypes)
    W:
        Matrix of archetypes (# features x # archetypes)
    H:
        Convex matrix of observation coefficients (# archetypes x # observations)    
        
    """
    
    AA_out = an.run_AA(A, k)            
    return(AA_out["C"], AA_out["H"])


def run_ACTION(
    data: Union[AnnData, np.ndarray],
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    thread_no: Optional[int] = 8,
    max_it: Optional[int] = 50,
    min_delta: Optional[float] = 0.01
) -> [dict]:
    """\
    Run ACTION decomposition [Mohammadi2018]_.

    Computes reduced ACTION decomposition.

    Parameters
    ----------
    data:
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    k_min:
        Min. # of archetypes to consider
    k_max:
        Max. # of archetypes to consider
    max_it, min_delta:
        Define stopping conditions of inner AA loop

    Returns
    -------
    C, H:
        A dictionary with trace of C and H matrices
    """
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        if 'ACTION_S_r' not in data.obsm.keys():
            print("ACTION_S_r is not in data.obsm. Please run reduce_kernel() first.")
            return ()
        X = data.obsm['ACTION_S_r'].T
    else:
        X = data.T

    print("Running ACTION ... ")
    ACTION_out = an.run_ACTION(X, k_min, k_max, thread_no, max_it, min_delta)
    print("Done.")
    
    return(ACTION_out["C"], ACTION_out["H"])
    
    

def prune_archetypes(
    ACTIONet_out: AnnData,
    C_trace: list,
    H_trace: list,
    min_specificity_z_threshold: Optional[float] = -1,
    copy: Optional[bool] = False    
) -> Optional[AnnData]:
    """\
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    ACTIONet_out:
        Current AnnData object storing the ACTIONet results    
    C_trace, H_trace:
        Output of run_ACTION()
    min_specificity_z_threshold:
        Controls level of prunning for non-specific archetypes (larger values remove more archetypes)..
    copy
        Determines whether a copy of ACTIONet_out is returned. 
    Returns
    -------
        None, if copy is False, and ACTIONet_out: AnnData, otherwise. 
        In either case, "C_stacked" and "H_stacked" are added to the ACTIONet_out.obsm.
    """

    print("Running archetype pruning ...")    
    prune_out = an.prune_archetypes(C_trace, H_trace, min_specificity_z_threshold)
    print("Done")
    
    ACTIONet_out.obsm['C_stacked'] = prune_out["C_stacked"]
    ACTIONet_out.obsm['H_stacked'] = prune_out["H_stacked"].T
    
    ACTIONet_out.uns["ACTIONet_pruning"] = {}
    reg = ACTIONet_out.uns["ACTIONet_pruning"]
    reg['selected_archs'] = prune_out["selected_archs"]    
    
    return ACTIONet_out if copy else None


def unify_archetypes(
    ACTIONet_out: AnnData,
    minPoints: Optional[int] = 10,
    minClusterSize: Optional[int] = 10,
    outlier_threshold: Optional[float] = 0, 
    copy: Optional[bool] = False       
) -> Optional[AnnData]:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    ACTIONet_out:
        Current AnnData object storing the ACTIONet results    
    int minPoints, int minClusterSize, double outlier_threshold:
        HDBSCAN parameters
    copy
        Determines whether a copy of ACTIONet_out is returned. 
    Returns
    -------
        None, if copy is False, and ACTIONet_out: AnnData, otherwise. 
        In either case, "C_pruned" and "H_pruned" are added to the ACTIONet_out.obsm.
    """
    if 'ACTION_S_r' not in ACTIONet_out.obsm.keys():
        print("ACTION_S_r is not in ACTIONet_out.obsm. Please run reduce_kernel() first.")
        return ACTIONet_out if copy else None
    else:
        S_r = ACTIONet_out.obsm["ACTION_S_r"].T
      
    if 'ACTIONet_connectivities' not in ACTIONet_out.obsp.keys():
        print("ACTIONet_connectivities is not in ACTIONet_out.obsp. Please run build_ACTIONet() first.")
        return ACTIONet_out if copy else None
    else:
        G = ACTIONet_out.obsp["ACTIONet_connectivities"]
        
    if 'H_stacked' not in ACTIONet_out.obsm.keys():
        print("H_stacked is not in ACTIONet_out.obsm. Please run prune_archetypes() first.")
        return ACTIONet_out if copy else None
    else:
        C_stacked = ACTIONet_out.obsm['C_stacked']
        H_stacked = ACTIONet_out.obsm['H_stacked'].T  
    

    
    unification_out = an.unify_archetypes(G, S_r, C_stacked, H_stacked, minPoints, minClusterSize, outlier_threshold)
    
    
    ACTIONet_out.obsm['C_unified'] = unification_out["C_unified"]
    ACTIONet_out.obsm['H_unified'] = unification_out["H_unified"].T
    
    ACTIONet_out.uns["ACTIONet_unification"] = {}
    reg = ACTIONet_out.uns["ACTIONet_unification"]
    reg['minPoints'] = minPoints
    reg['minClusterSize'] = minClusterSize
    reg['outlier_threshold'] = outlier_threshold
    reg['archetype_groupd'] = unification_out["archetype_groups"]

    groups = unification_out["sample_assignments"]
    ACTIONet_out.obs["assigned_archetype"] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )    
        
    
    return ACTIONet_out if copy else None
