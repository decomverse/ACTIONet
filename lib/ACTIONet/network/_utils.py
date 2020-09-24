from typing import Optional

from anndata import AnnData

import _ACTIONet as _an

def compute_archetype_core_centrality(
    adata: AnnData,
    key: Optional[str] = 'ACTIONet',    
    copy: Optional[bool] = False    
) -> AnnData:
    """\
    Computes node centrality scores
    
    Uses graph core-ness to compute node centralities
    
    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    key
        `adata.obsp` key that stores the ACTIONet connectivities
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, ACE: AnnData, if copy is True.        
        "node_centrality" is to ACE.obs.        
    """    
    if 'ACTIONet' not in adata.obsp.keys():
        raise ValueError(
            'Did not find adata.obsp[\'ACTIONet\']. '
            'Please run nt.build_network() first.'
        )
    if 'ACTION' not in adata.obs.keys():
        raise ValueError(
            'Did not find adata.obs[\'ACTION\']. '
            'Please run pp.unify_archetypes() first.'
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp['ACTIONet']
    assignments = adata.obs['ACTION']

    scores = _an.compute_archetype_core_centrality(G, assignments)
    adata.obs['ACTIONet_centrality'] = scores
    
    return adata if copy else None
