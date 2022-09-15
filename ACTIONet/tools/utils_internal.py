from typing import Optional, Union

import pandas as pd
from anndata import AnnData

# def __get_attr_or_split_idx(
#         adata: AnnData,
#         attr: Union[str, list, pd.Series],
#         groups_use: Union[str, list, None] = None,
#         return_vec: Optional[bool] = False,
#         d: Optional[int] = 0,
#         ) -> Union[list, dict]:
#     if d not in [0, 1]:
#         raise ValueError("d must be dim (0 or 1) of adata")
#
#     if isinstance(attr, str) or len(attr) == 1:
#         if d == 0:
#             split_vec = adata.obs[attr]
#         else:
#             split_vec = adata.var[attr]
#
#     else:
#         if len(attr) != adata.shape[d]:
#             raise ValueError("len(attr) does not match .shape[{dim:d}] of adata".format(dim=d))
#         split_vec = attr
#
#     if split_vec is None:
#         raise ValueError("Invalid split condition")
#     else:
#         if isinstance(split_vec, pd.Series):
#             split_vec = split_vec.tolist()
#
#     idx = list(range(0, adata.shape[d]))
#     idx_out = {}
#
#     if groups_use is not None:
#         idx = [i for i, e in enumerate(split_vec) if e in groups_use]
#         split_vec = [e for i, e in enumerate(split_vec) if e in groups_use]
#
#         if len(split_vec) == 0:
#             raise ValueError("Invalid split condition")
#
#         for g in set(groups_use):
#             idx_out[g] = [idx.index(i) for i, e in enumerate(split_vec) if e == g]
#
#     else:
#         for g in set(split_vec):
#             idx_out[g] = [idx.index(i) for i, e in enumerate(split_vec) if e == g]
#
#     return split_vec if return_vec else idx_out


def __get_feature_vec(
    adata: AnnData,
    features_use: Union[str, list, pd.Series, None] = None,
) -> Union[list, dict]:
    if features_use is None:
        features_use = adata.var_names
    else:
        # features_use = __get_attr_or_split_idx(
        #         adata=adata,
        #         attr=features_use,
        #         return_vec=True,
        #         d=1
        #         )
        features_use = get_data_or_split(
            adata=adata, attr=features_use, groups_use=None, to_return="data", d=1
        )
    return features_use
