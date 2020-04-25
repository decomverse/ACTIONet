import _ACTIONet as an
from typing import Optional
from anndata import AnnData

    
def build_ACTIONet(
    ACTIONet_out: AnnData,
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
    ACTIONet_out:
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
        None, if copy is False, and ACTIONet_out: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACTIONet_out.obsp.
    """

    if 'H_stacked' not in ACTIONet_out.obsm.keys():
        print("H_stacked is not in ACTIONet_out.obsm. Please run prune_archetypes() first.")
        return ACTIONet_out if copy else None
    else:
        H_stacked = ACTIONet_out.obsm["H_stacked"].T
  
    

    print("Building ACTIONet graph ...")
    G = an.build_ACTIONet(H_stacked, density, thread_no, mutual_edges_only)
    print("Done.")
    
    ACTIONet_out.uns["ACTIONet_build"] = {}
    neighbors_dict = ACTIONet_out.uns["ACTIONet_build"]
    neighbors_dict['density'] = density
    neighbors_dict['mutual_edges_only'] = mutual_edges_only

    ACTIONet_out.obsp["ACTIONet_connectivities"] = G

    return ACTIONet_out if copy else None


def layout_ACTIONet(
    ACTIONet_out: AnnData,
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
    ACTIONet_out:
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
        None, if copy is False, ACTIONet_out: AnnData, if copy is True.        
        "ACTIONet2D" and "ACTIONet3D" are added to the ACTIONet_out.obsm, and "denovo_colors" to ACTIONet_out.uns.
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
      
    
    vis_out = an.layout_ACTIONet(G, S_r, compactness_level, n_epochs, thread_no)
        
    ACTIONet_out.obsm["X_ACTIONet2D"] = vis_out["coordinates"]
    ACTIONet_out.obsm["X_ACTIONet3D"] = vis_out["coordinates_3D"]
    ACTIONet_out.uns["denovo_colors"] = vis_out["colors"]
    
    return ACTIONet_out if copy else None


def compute_archetype_core_centrality(
    ACTIONet_out: AnnData,
    key: Optional[str] = "ACTIONet_connectivities",    
    copy: Optional[bool] = False    
) -> Optional[AnnData]:
    """\
    Computes node centrality scores
    
    Uses graph core-ness to compute node centralities
    
    Parameters
    ----------
    ACTIONet_out:
        AnnData object storing the ACTIONet results        
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, ACTIONet_out: AnnData, if copy is True.        
        "node_centrality" is to ACTIONet_out.obs.        
    """    
    
    if 'ACTIONet_connectivities' not in ACTIONet_out.obsp.keys():
        print("ACTIONet_connectivities is not in ACTIONet_out.obsp. Please run build_ACTIONet() first.")
        return ACTIONet_out if copy else None
    else:
        G = ACTIONet_out.obsp["ACTIONet_connectivities"]

    if 'assigned_archetype' not in ACTIONet_out.obs.keys():
        print("assigned_archetype is not in ACTIONet_out.obs. Please run unify_archetypes() first.")
        return ACTIONet_out if copy else None
    else:
        assigned_archetype = ACTIONet_out.obs['assigned_archetype']

        
    scores = an.compute_archetype_core_centrality(G, assigned_archetype)
    
    ACTIONet_out.obs["node_centrality"] = scores
    return ACTIONet_out if copy else None
