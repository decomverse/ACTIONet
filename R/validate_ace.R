## Validates an ACTIONetExperiment (ACE) object
#' @importFrom BiocGenerics NCOL NROW
setValidity2("ACTIONetExperiment", function(object) {
    NR <- NROW(object)
    NC <- NCOL(object)
    msg <- NULL

    # 2D
    value = object@rowNets
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NR) | (NCOL(value[[i]]) != NR)) {
            msg <- c(msg, "'nrow(rowNets[[...]])' and 'ncol(rowNets[[...]])' should be equal to the number of rows of ace object.")
        }
    }

    value = object@colNets
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NC) | (NCOL(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colNets[[...]])' and 'ncol(colNets[[...]])' should be equal to the number of columns of ace object.")
        }
    }

    value = object@rowMaps
    for (i in seq_along(value)) {
        value[[i]] = .validate_MapType(value[[i]])
        if ((NROW(value[[i]]) != NR)) {
            msg <- c(msg, "'nrow(rowMaps[[..]])' should be equal to the number of rows of ace object..")
        }

        if (any(rownames(value[[i]]) != rownames(object))) {
            msg <- c(msg, "'rownames(rowMaps[[..]])' must match the rownames of ace object.")
        }
    }


    value = object@colMaps
    for (i in seq_along(value)) {
        value[[i]] = .validate_MapType(value[[i]])
        if ((NROW(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colMaps[[..]])' should be equal to the number of columns of ace object..")
        }

        if (any(rownames(value[[i]]) != colnames(object))) {
            msg <- c(msg, "'rownames(colMaps[[..]])' must match the colnames of ace object.")
        }
    }


    if (length(msg)) {
        msg
    } else TRUE
})

.validate_MapType <- function(SE) {
    value = S4Vectors::metadata(SE)$type

    if (!(value %in% c("generic", "reduction", "embedding", "internal")) | is.null(value)) {
        err = sprintf("MapType must be 'generic', 'reduction', 'embedding', or 'internal'. Setting to 'generic'.\n")
        warning(err)
        S4Vectors::metadata(SE)$type = "generic"
        return(SE)
    }
    return(SE)
}

.validate_names <- function(value, valid_names = NULL) {
    par_func = as.character(sys.call(-1)[1])
    if (any(names(value) == "")) {
        err = sprintf("Values passed to '%s' cannot be unnamed.\n", par_func)
        stop(err)
    }

    if (any(duplicated(names(value)))) {
        err = sprintf("Values passed to '%s' have duplicate names.\n", par_func)
        stop(err)
    }

    if (!is.null(valid_names)) {
        not_in_object = setdiff(names(value), valid_names)
        if (length(not_in_object) > 0) {
            err = sprintf("No element named '%s'.\n", not_in_object)
            stop(err)
        }
    }

    return
}
