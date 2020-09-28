from typing import Optional

from anndata import AnnData

import _ACTIONet as _an

def build_network(
    adata: AnnData,
    density: Optional[float] = 1.0,
    n_threads: Optional[int] = 0,
    mutual_edges_only: Optional[bool] = True,
    copy: Optional[bool] = False        
) -> AnnData:
    """\
    Build ACTIIONet 
    
    Computes and returns the ACTIONet graph
    Parameters
    ----------
    adata
        Current AnnData object storing the ACTIONet results    
    density
        Controls the overall density of constructed network.
        Larger values results in more retained edges.
    n_threads
        Number of parallel threads used for identifying nearest-neighbors.
        Defaults to available threads on the machine.
    mutual_edges_only
        Whether to return only edges that there is a bi-directional/mutual relationship
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsp['ACTIONet']`
        `.uns['ACTIONet']['params']`
    """
    if 'ACTION_H_stacked' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'ACTION_H_stacked\']. '
            'Please run pp.prune_archetypes() first.'
        )
    adata = adata.copy() if copy else adata
    H_stacked = adata.obsm["ACTION_H_stacked"].T
    G = _an.build_ACTIONet(H_stacked, density, n_threads, mutual_edges_only)
    
    adata.uns.setdefault('ACTIONet', {}).update({'params': {
        'density': density,
        'mutual_edges_only': mutual_edges_only,
    }})
    adata.obsp['ACTIONet'] = G

    return adata if copy else None

