from typing import Union

import pandas as pd
from anndata import AnnData


def __get_feature_vec(
    adata: AnnData,
    features_use: Union[str, list, pd.Series, None] = None,
) -> Union[list, dict]:
    return adata.var_names
