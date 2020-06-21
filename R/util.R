#' @export
get_ace_split_IDX <- function(ace, attr, groups_use = NULL, return_split_vec = FALSE) {
    if (length(attr) == 1) {
        split_vec = colData(ace)[[attr]]
    } else {
        split_vec = attr
    }
    
    col_idx = 1:dim(ace)[2]
    
    if (!is.null(groups_use)) {
        sub_idx = which(split_vec %in% groups_use)
        split_vec = split_vec[sub_idx]
        col_idx = col_idx[sub_idx]
    }
    
    if (is.null(split_vec)) 
        stop(sprintf("Invalid split conditions.\n"))
    
    split_vec = droplevels(factor(split_vec))
    IDX = split(col_idx, split_vec)
    
    if (return_split_vec) 
        return(split_vec) else return(IDX)
}

check_if_ace <- function(sce_like) {
    if (class(sce_like) %in% c("ACTIONetExperiment", "SummarizedExperiment", "SingleCellExperiment")) {
        if (class(sce_like) != "ACTIONetExperiment") {
            ace = as(sce_like, "ACTIONetExperiment")
            msg = sprintf("Converting to ACTIONetExperiment class.\n")
            message(msg)
            return(ace)
        } else {
            return(sce_like)
        }
    } else {
        err = sprintf("Input must type ACTIONetExperiment, SingleCellExperiment, or SummarizedExperiment.\n")
        stop(err)
    }
}
