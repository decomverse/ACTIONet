import _ACTIONet as an
from typing import Optional
from anndata import AnnData

    
def build_ACTIONet(
    ACE: AnnData,
    density: Optional[float] = 1,
    thread_no: Optional[int] = 8, 
    mutual_edges_only: Optional[bool] = True,
    copy: Optional[bool] = False        
) -> Optional[AnnData]:
    """\
    Build ACTIIONet 
    
    Computes and returns the ACTIONet graph
    Parameters
    ----------
    ACE:
        Current AnnData object storing the ACTIONet results    
    density:
        Controls the overall density of constructed network. Larger values results in more retained edges.
    thread_no:
        Number of parallel threads used for identifying nearest-neighbors
    mutual_edges_only:
        Whether to return only edges that there is a bi-directional/mutual relationship
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, and ACE: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACE.obsp.
    """

    if 'H_stacked' not in ACE.obsm.keys():
        print("H_stacked is not in ACE.obsm. Please run prune_archetypes() first.")
        return ACE if copy else None
    else:
        H_stacked = ACE.obsm["H_stacked"].T
  
    

    print("Building ACTIONet graph ...")
    G = an.build_ACTIONet(H_stacked, density, thread_no, mutual_edges_only)
    print("Done.")
    
    ACE.uns["ACTIONet_build"] = {}
    neighbors_dict = ACE.uns["ACTIONet_build"]
    neighbors_dict['density'] = density
    neighbors_dict['mutual_edges_only'] = mutual_edges_only

    ACE.obsp["ACTIONet_connectivities"] = G

    return ACE if copy else None


def layout_ACTIONet(
    ACE: AnnData,
    compactness_level: Optional[int] = 50,
    n_epochs: Optional[int] = 500,
    thread_no: Optional[int] = 8,
    copy: Optional[bool] = False
) -> Optional[AnnData]:
    """\
    Network layout
    
    Embedded the graph into 2D/3D space
    
    Parameters
    ----------
    ACE:
        AnnData object storing the ACTIONet results
    compactness_level:
        Between 0-100
    n_epochs:
        Number of SGD epochs
    thread_no:
        Number of parallel threads
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, ACE: AnnData, if copy is True.        
        "ACTIONet2D" and "ACTIONet3D" are added to the ACE.obsm, and "denovo_colors" to ACE.uns.
    """
    if 'ACTION_S_r' not in ACE.obsm.keys():
        print("ACTION_S_r is not in ACE.obsm. Please run reduce_kernel() first.")
        return ACE if copy else None
    else:
        S_r = ACE.obsm["ACTION_S_r"].T
    
    if 'ACTIONet_connectivities' not in ACE.obsp.keys():
        print("ACTIONet_connectivities is not in ACE.obsp. Please run build_ACTIONet() first.")
        return ACE if copy else None
    else:
        G = ACE.obsp["ACTIONet_connectivities"]
      
    
    vis_out = an.layout_ACTIONet(G, S_r, compactness_level, n_epochs, thread_no)
        
    ACE.obsm["X_ACTIONet2D"] = vis_out["coordinates"]
    ACE.obsm["X_ACTIONet3D"] = vis_out["coordinates_3D"]
    ACE.uns["denovo_colors"] = vis_out["colors"]
    
    return ACE if copy else None


def compute_archetype_core_centrality(
    ACE: AnnData,
    key: Optional[str] = "ACTIONet_connectivities",    
    copy: Optional[bool] = False    
) -> Optional[AnnData]:
    """\
    Computes node centrality scores
    
    Uses graph core-ness to compute node centralities
    
    Parameters
    ----------
    ACE:
        AnnData object storing the ACTIONet results        
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, ACE: AnnData, if copy is True.        
        "node_centrality" is to ACE.obs.        
    """    
    
    if 'ACTIONet_connectivities' not in ACE.obsp.keys():
        print("ACTIONet_connectivities is not in ACE.obsp. Please run build_ACTIONet() first.")
        return ACE if copy else None
    else:
        G = ACE.obsp["ACTIONet_connectivities"]

    if 'assigned_archetypes' not in ACE.obs.keys():
        print("assigned_archetypes is not in ACE.obs. Please run unify_archetypes() first.")
        return ACE if copy else None
    else:
        assigned_archetypes = ACE.obs['assigned_archetypes']

        
    scores = an.compute_archetype_core_centrality(G, assigned_archetypes)
    
    ACE.obs["node_centrality"] = scores
    return ACE if copy else None
