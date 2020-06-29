import _ACTIONet as an
import pandas as pd
from scipy.sparse import issparse
from typing import Optional
from anndata import AnnData

    
def compute_archetype_feature_specificity(
    ACE: AnnData,
    archetype_key: Optional[str] = "H_unified",
    copy: Optional[bool] = False        
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of archetypes
    
    Uses Archetype footprints to estimate markers (soft clustering)
    Parameters
    ----------
    ACE:
        Current AnnData object storing the ACTIONet results   
    archetype_key:
        Key in ACE.obsm that holds the archetype footprints (default = "H_unified")
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, and ACE: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACE.obsp.
    """

    S = ACE.X.T
    
    if archetype_key not in ACE.obsm.keys():
        print(archetype_key + " is not in ACE.obsm.")
        return ACE if copy else None
    else:
        H_stacked = ACE.obsm[archetype_key].T
  
    
    if (issparse(S)):    
        print("Running archetype feature specificity assessment in sparse mode ...")    
        spec_out = an.compute_archetype_feature_specificity (S, H_stacked)
    else:
        print("Running archetype feature specificity assessment in dense mode ...")    
        spec_out = an.compute_archetype_feature_specificity_full(S, H_stacked)
        
    ACE.varm[archetype_key+"_profile"] = spec_out["archetypes"]
    ACE.varm[archetype_key+"_upper_significance"] = spec_out["upper_significance"]
    ACE.varm[archetype_key+"_lower_significance"] = spec_out["lower_significance"]
    
    print("Done.")    
        
    return ACE if copy else None


def compute_cluster_feature_specificity(
    ACE: AnnData,
    cluster_key: Optional[str] = "leiden",
    copy: Optional[bool] = False        
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of clusters
    
    Uses cluster membership vector to estimate markers (disjoint clustering)
    Parameters
    ----------
    ACE:
        Current AnnData object storing the ACTIONet results   
    cluster_key:
        Key in ACE.obs that holds the clustering variable (default = "leiden")
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, and ACE: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACE.obsp.
    """

    S = ACE.X.T
    if cluster_key not in ACE.obs.keys():
        print(cluster_key + " is not in ACE.obsm.")
        return ACE if copy else None
    else:
        clusters = ACE.obs[cluster_key]
  
    if isinstance(clusters, pd.Series):
        clusters = pd.factorize(clusters)[0]
    
    if (issparse(S)):    
        print("Running cluster feature specificity assessment in sparse mode ...")    
        spec_out = an.compute_cluster_feature_specificity (S, clusters)
    else:
        print("Running cluster feature specificity assessment in dense mode ...")    
        spec_out = an.compute_cluster_feature_specificity_full(S, clusters)
        
    ACE.varm[cluster_key+"_profile"] = spec_out["archetypes"]
    ACE.varm[cluster_key+"_upper_significance"] = spec_out["upper_significance"]
    ACE.varm[cluster_key+"_lower_significance"] = spec_out["lower_significance"]
    
    print("Done.")    
        
    return ACE if copy else None
