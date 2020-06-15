#' @export
get_ace_split_IDX <- function(ace, attr, groups_use = NULL){
  if( length(attr) ==  1) {
    split_vec = colData(ace)[[attr]]
  }
  else {
    split_vec = attr
  }

  col_idx = 1:dim(ace)[2]

  if(!is.null(groups_use)){
    sub_idx = which(split_vec %in% groups_use)
    split_vec = split_vec[sub_idx]
    col_idx = col_idx[sub_idx]
  }

  if(is.null(split_vec))
    stop(sprintf("Invalid split conditions.\n"))

  split_vec = droplevels(factor(split_vec))
  IDX = split(col_idx, split_vec)

  return(IDX)
}
