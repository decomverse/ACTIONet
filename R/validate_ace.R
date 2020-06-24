#' Validates an ACTIONetExperiment (ACE) object
#'
#' @param object ACTIONetExperiment object
#'
#' @importFrom BiocGenerics NCOL NROW
setValidity2("ACTIONetExperiment", function(object) {
    NR <- NROW(object)
    NC <- NCOL(object)
    msg <- NULL

    # 2D
    value = rowNets(object)
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NR) | (NCOL(value[[i]]) != NR)) {
            msg <- c(msg, "'nrow(rowNets[[...]])' and 'ncol(rowNets[[...]])' should be equal to the number of rows of ace object.")
        }
    }

    value = colNets(object)
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NC) | (NCOL(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colNets[[...]])' and 'ncol(colNets[[...]])' should be equal to the number of columns of ace object.")
        }
    }

    value = rowMaps(object)
    for (i in seq_along(value)) {
        if ((NCOL(value[[i]]) != NR)) {
            msg <- c(msg, "'ncol(rowMaps[[..]])' should be equal to the number of rows of ace object..")
        }

        if (any(colnames(value[[i]]) != rownames(object))) {
            msg <- c(msg, "'colnames(rowMaps[[..]])' must match the rownames of ace object.")
        }
    }


    value = colMaps(object)
    for (i in seq_along(value)) {
        if ((NCOL(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colMaps[[..]])' should be equal to the number of columns of ace object..")
        }

        if (any(colnames(value[[i]]) != colnames(object))) {
            msg <- c(msg, "'colnames(colMaps[[..]])' must match the colnames of ace object.")
        }
    }

    if (length(msg)) {
        msg
    } else TRUE
})
