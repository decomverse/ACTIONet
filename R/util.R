.get_ace_split_IDX <- function(ace, attr = NULL){
  if( length(attr) ==  1) {
    split_vec = colData(sce)[[attr]]
  }
  else {
    split_vec = attr
  }
  split_vec = droplevels(factor(split_vec))
  IDX = split(1:dim(sce)[2], split_vec)

  return(IDX)
}
