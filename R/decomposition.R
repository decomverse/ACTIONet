#' @export
decomp <- function(X, k, method, W0 = NULL, H0 = NULL, params) {
    if (method == "simplex_regression") {
        W <- as.matrix(W0)
        B <- as.matrix(X)
        H <- run_simplex_regression(A = W, B = B)
        extra <- list()
    } else if (method == "SVD") {
        if (is.matrix(X)) {
            reduction.out <- IRLB_SVD_full(
                A = X, dim = k,
                iter = params$max_iter, seed = params$seed
            )
        } else {
            reduction.out <- IRLB_SVD(
                A = X, dim = k,
                iter = params$max_iter, seed = params$seed
            )
        }
        W <- reduction.out$u
        H <- as.matrix(Diagonal(length(reduction.out$d), reduction.out$d) %*% t(reduction.out$v))
        extra <- list(d = reduction.out$d)
    } else if (method == "PCA") {
        SVD.out <- decomp(X, k, "SVD", params = params)
        u <- as.matrix(SVD.out$W)
        v <- Matrix::t(as.matrix(Diagonal(length(SVD.out$extra$d), 1 / SVD.out$extra$d) %*% SVD.out$H))
        d <- as.matrix(SVD.out$extra$d)
        if (is.matrix(X)) {
            reduction.out <- SVD2PCA_full(S = X, u = u, d = d, v = v)
        } else {
            reduction.out <- SVD2PCA(X = X, u = u, d = d, v = v)
        }
        W <- reduction.out$x
        H <- reduction.out$rotation
        extra <- list(sdev = reduction.out$sdev)
        extra$d <- extra$sdev * sqrt(nrow(W) - 1)
    } else if (method == "ACTION_reduction") {
        if (is.matrix(X)) {
            reduction.out <- reduce_kernel_full(
                S = X, reduced_dim = k,
                iter = params$max_iter, seed = params$seed, SVD_algorithm = params$SVD_algorithm
            )
        } else {
            reduction.out <- reduce_kernel(
                S = X, reduced_dim = k,
                iter = params$max_iter, seed = params$seed, SVD_algorithm = params$SVD_algorithm
            )
        }
        W <- reduction.out$V
        H <- reduction.out$S_r
        extra <- list(A = reduction.out$A, B = reduction.out$B, reduction.out$sigma)
    } else if (method == "SPA") {
        out <- run_SPA(X, k)
        W <- X[, out$selected_columns]
        H <- run_simplex_regression(W, X)

        extra <- list(selected_columns = out$selected_columns, norms = out$norms)
    } else if (method == "AA") {
        if (is.null(W0)) {
            W0 <- matrix(0, nrow(X), k)
        }
        out <- run_AA(as.matrix(X), W0 = as.matrix(W0))
        W <- out$W
        H <- out$H
        extra <- list(C = out$C)
    } else if (method == "ACTION_decomposition") {
        out <- run_ACTION(S_r = X, k_min = k, k_max = k, thread_no = params$thread_no, max_it = params$max_it, min_delta = 1e-300)
        H <- out$H[[k]]
        C <- out$C[[k]]
        W <- X %*% C
        extra <- list(C = C)
    } else if (method == "ACTION_decomposition_MR") {
        ACTION.out <- run_ACTION(
            S_r = X,
            k_min = params$k_min,
            k_max = params$k_max,
            thread_no = params$thread_no,
            max_it = params$max_iter,
            min_delta = 1e-300
        )

        # Prune nonspecific and/or unreliable archetypes
        pruning.out <- prune_archetypes(
            C_trace = ACTION.out$C,
            H_trace = ACTION.out$H,
            min_specificity_z_thresh = params$min_specificity_z_thresh,
            min_cells = params$min_cells_per_arch
        )

        # Identiy equivalent classes of archetypes and group them together
        unification.out <- unify_archetypes(
            S_r = X,
            C_stacked = pruning.out$C_stacked,
            H_stacked = pruning.out$H_stacked,
            violation_threshold = params$unification_violation_threshold,
            thread_no = params$thread_no
        )
        H <- unification.out$H_unified
        W <- X %*% unification.out$C_unified
        extra <- list(H = ACTION.out$H, C = ACTION.out$C, H_stacked = pruning.out$H_stacked, C_stacked = pruning.out$C_stacked, H_unified = unification.out$H_unified, C_unified = unification.out$C_unified, assigned_archetype = unification.out$assigned_archetype)
    } else if (method == "ACTION_decomposition_MR_with_batch_correction") {
        ACTION.out <- run_ACTION_with_batch_correction(
            S_r = X,
            batch = params$batch,
            k_min = params$k_min,
            k_max = params$k_max,
            thread_no = params$thread_no,
            max_it = params$max_iter,
            min_delta = 1e-300
        )

        # Prune nonspecific and/or unreliable archetypes
        pruning.out <- prune_archetypes(
            C_trace = ACTION.out$C,
            H_trace = ACTION.out$H,
            min_specificity_z_thresh = params$min_specificity_z_thresh,
            min_cells = params$min_cells_per_arch
        )

        # Identiy equivalent classes of archetypes and group them together
        unification.out <- unify_archetypes(
            S_r = X,
            C_stacked = pruning.out$C_stacked,
            H_stacked = pruning.out$H_stacked,
            violation_threshold = params$unification_violation_threshold,
            thread_no = params$thread_no
        )
        H <- unification.out$H_unified
        W <- X %*% unification.out$C_unified
        extra <- list(H = ACTION.out$H, C = ACTION.out$C, H_stacked = pruning.out$H_stacked, C_stacked = pruning.out$C_stacked, H_unified = unification.out$H_unified, C_unified = unification.out$C_unified, assigned_archetype = unification.out$assigned_archetype)
    }

    out <- list(W = W, H = H, extra = extra)
    return(out)
}