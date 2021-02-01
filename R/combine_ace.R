
#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("cbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    # args <- lapply(args, updateObject)
    args <- .clear_ace_slots_bind(args)
    out <- do.call(callNextMethod, args)
    return(out)
})

#' @importFrom SummarizedExperiment rbind cbind
#' @export
setMethod("rbind", "ACTIONetExperiment", function(..., deparse.level = 1) {
    args <- list(...)
    # args <- lapply(args, updateObject)
    args <- .clear_ace_slots_bind(args)
    out <- do.call(callNextMethod, args)
    return(out)
})

.clear_ace_slots_bind <- function(args) {
    used_slots = c()
    AN_slots = c("colMaps", "rowMaps", "colNets", "rowNets")
    used_slots = sapply(args, function(x) {
        sl = c(length(colMaps(x)) > 0, length(rowMaps(x)) > 0, length(colNets(x)) > 
            0, length(rowNets(x)) > 0)
        return(AN_slots[sl])
    })
    
    used_slots = Reduce(union, c(used_slots))
    if (length(used_slots) > 0) {
        par_func = as.character(sys.call(-1)[1])
        w = paste(sprintf("In %s: ", par_func), "Non-concatable slot <(", used_slots, 
            sprintf(")> will not be preserved.\n"), sep = "")
        warning(w, call. = FALSE)
    }
    
    nc_rep = S4Vectors::SimpleList()
    args = lapply(args, function(a) {
        BiocGenerics:::replaceSlots(a, rowNets = nc_rep, colNets = nc_rep, rowMaps = nc_rep, 
            colMaps = nc_rep, check = FALSE)
    })
    return(args)
}
