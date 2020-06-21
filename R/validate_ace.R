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
            msg <- c(msg, "'nrow(rowNets[[...]])' and 'ncol(rowNets[[...]])' should be equal to the number of rows")
        }
    }

    value = colNets(object)
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NC) | (NCOL(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colNets[[...]])' and 'ncol(colNets[[...]])' should be equal to the number of columns")
        }
    }

    value = rowMaps(object)
    for (i in seq_along(value)) {
        if ((NROW(value[[i]]) != NR)) {
            msg <- c(msg, "'nrow(rowMaps[[..]])' should be equal to the number of rows.")
        }

        if ((rownames(value[[i]]) != rownames(object))) {
            msg <- c(msg, "'rownames(rowMaps[[..]])' must match the rownames of ace object.")
        }
    }


    value = colMaps(object)
    for (i in seq_along(value)) {
        if ((NCOL(value[[i]]) != NC)) {
            msg <- c(msg, "'nrow(colMaps[[..]])' should be equal to the number of columns.")
        }

        if ((colnames(value[[i]]) != colnames(object))) {
            msg <- c(msg, "'colnames(colMaps[[..]])' must match the colnames of ace object.")
        }
    }

    if (length(msg)) {
        msg
    } else TRUE
})
