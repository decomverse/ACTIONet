from typing import Optional, Union
from anndata import AnnData
import numpy as np
from scipy import sparse

import _ACTIONet as _an


def build_network(
        data: Union[AnnData, np.ndarray, sparse.spmatrix],
        net_name_out: Optional[str] = "ACTIONet",
        density: Optional[float] = 1.0,
        thread_no: Optional[int] = 0,
        mutual_edges_only: Optional[bool] = True,
        distance_metric:str = "jsd",
        nn_approach:str = "k*nn",
        k: Optional[int] = None,
        M: Optional[float] = 16,
        ef_construction:Optional[float] = 200,
        ef:Optional[float] = 50,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False,
        
        
) -> Union[AnnData, sparse.spmatrix, None]:

    """[Summary]

    :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
    :type [ParamName]: [ParamType](, optional)
    ...
    :raises [ErrorType]: [ErrorDescription]
    ...
    :return: [ReturnDescription]
    :rtype: [ReturnType]
    """

    data_is_AnnData = isinstance(data, AnnData)

    if data_is_AnnData:
        adata = data.copy() if copy else data
        if "H_stacked" in adata.obsm:
            H_stacked = adata.obsm["H_stacked"]
        else:
            raise Exception("missing H_stacked key in adata.obsm of AnnData")
    else:
        H_stacked=data.X 
    H_stacked = H_stacked.T.astype(dtype=np.float64)
    if sparse.issparse(H_stacked):
        H_stacked = H_stacked.toarray()
        if distance_metric=="jsd":
            H_stacked=H_stacked/np.sum(H_stacked,axis=0)
        G = _an.build_ACTIONet(H_stacked=H_stacked,
                               density=density,
                               thread_no=thread_no,
                               mutual_edges_only=mutual_edges_only,
                               distance_metric=distance_metric,
                               nn_approach=nn_approach,
                               k=k,
                               M=M,
                               ef_construction=ef_construction,
                               ef=ef)

    if return_raw or not data_is_AnnData:
        return G
    elif data_is_AnnData:
        adata.obsp[net_name_out] = G
        return adata if copy else None

