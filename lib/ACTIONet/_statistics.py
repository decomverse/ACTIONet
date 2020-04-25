import _ACTIONet as an
import pandas as pd
from scipy.sparse import issparse
from typing import Optional
from anndata import AnnData

    
def compute_archetype_feature_specificity(
    ACTIONet_out: AnnData,
    archetype_key: Optional[str] = "H_unified",
    copy: Optional[bool] = False        
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of archetypes
    
    Uses Archetype footprints to estimate markers (soft clustering)
    Parameters
    ----------
    ACTIONet_out:
        Current AnnData object storing the ACTIONet results   
    archetype_key:
        Key in ACTIONet_out.obsm that holds the archetype footprints (default = "H_unified")
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, and ACTIONet_out: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACTIONet_out.obsp.
    """

    S = ACTIONet_out.X.T
    
    if archetype_key not in ACTIONet_out.obsm.keys():
        print(archetype_key + " is not in ACTIONet_out.obsm.")
        return ACTIONet_out if copy else None
    else:
        H_stacked = ACTIONet_out.obsm[archetype_key].T
  
    
    if (issparse(S)):    
        print("Running archetype feature specificity assessment in sparse mode ...")    
        spec_out = an.compute_archetype_feature_specificity_sparse (S, H_stacked)
    else:
        print("Running archetype feature specificity assessment in dense mode ...")    
        spec_out = an.compute_archetype_feature_specificity(S, H_stacked)
        
    ACTIONet_out.varm[archetype_key+"_profile"] = spec_out["archetypes"]
    ACTIONet_out.varm[archetype_key+"_upper_significance"] = spec_out["upper_significance"]
    ACTIONet_out.varm[archetype_key+"_lower_significance"] = spec_out["lower_significance"]
    
    print("Done.")    
        
    return ACTIONet_out if copy else None


def compute_cluster_feature_specificity(
    ACTIONet_out: AnnData,
    cluster_key: Optional[str] = "leiden",
    copy: Optional[bool] = False        
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of clusters
    
    Uses cluster membership vector to estimate markers (disjoint clustering)
    Parameters
    ----------
    ACTIONet_out:
        Current AnnData object storing the ACTIONet results   
    cluster_key:
        Key in ACTIONet_out.obs that holds the clustering variable (default = "leiden")
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, and ACTIONet_out: AnnData, otherwise. 
        In either case, "ACTIONet_connectivities" is added to ACTIONet_out.obsp.
    """

    S = ACTIONet_out.X.T
    if cluster_key not in ACTIONet_out.obs.keys():
        print(cluster_key + " is not in ACTIONet_out.obsm.")
        return ACTIONet_out if copy else None
    else:
        clusters = ACTIONet_out.obs[cluster_key]
  
    if isinstance(clusters, pd.Series):
        clusters = pd.factorize(clusters)[0]
    
    if (issparse(S)):    
        print("Running cluster feature specificity assessment in sparse mode ...")    
        spec_out = an.compute_cluster_feature_specificity_sparse (S, clusters)
    else:
        print("Running cluster feature specificity assessment in dense mode ...")    
        spec_out = an.compute_cluster_feature_specificity(S, clusters)
        
    ACTIONet_out.varm[cluster_key+"_profile"] = spec_out["archetypes"]
    ACTIONet_out.varm[cluster_key+"_upper_significance"] = spec_out["upper_significance"]
    ACTIONet_out.varm[cluster_key+"_lower_significance"] = spec_out["lower_significance"]
    
    print("Done.")    
        
    return ACTIONet_out if copy else None
