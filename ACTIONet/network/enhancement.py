import numpy as np
import scipy.sparse as sp


def dominateset(aff_matrix: np.ndarray, NR_OF_KNN: int) -> np.ndarray:
    """
    Compute the dominant nearest neighbors matrix.

    Args:
        aff_matrix: The affinity matrix.
        NR_OF_KNN: The number of nearest neighbors to consider.

    Returns:
        The dominant nearest neighbors matrix.
    """
    A = np.sort(aff_matrix, axis=1)[:, ::-1]
    res = A[:, :NR_OF_KNN]
    inds = np.repeat(np.arange(len(aff_matrix))[:, np.newaxis], NR_OF_KNN, axis=1)
    loc = np.argsort(-aff_matrix, axis=1)[:, :NR_OF_KNN]
    PNN_matrix1 = np.zeros_like(aff_matrix)
    PNN_matrix1[inds.flatten(), loc.flatten()] = res.flatten()
    PNN_matrix = (PNN_matrix1 + PNN_matrix1.T) / 2
    return PNN_matrix


def TransitionFields(W: np.ndarray) -> np.ndarray:
    """
    Compute the transition fields.

    Args:
        W: The input matrix.

    Returns:
        The transition fields matrix.
    """
    zeroindex = np.where(np.sum(W, axis=1) == 0)[0]
    W = W * len(W)
    W = dn(W, "ave")
    w = np.sqrt(np.sum(np.abs(W), axis=1) + np.finfo(float).eps)
    W = W / np.repeat(w[:, np.newaxis], len(W), axis=1)
    W = W @ W.T
    Wnew = W
    Wnew[zeroindex, :] = 0
    Wnew[:, zeroindex] = 0
    return Wnew


def dn(w: np.ndarray, type: str) -> np.ndarray:
    """
    Compute the normalized matrix.

    Args:
        w: The input matrix.
        type: The type of normalization to apply.

    Returns:
        The normalized matrix.
    """
    w = w * len(w)
    D = np.sum(np.abs(w), axis=1) + np.finfo(float).eps
    if type == "ave":
        D = 1 / D
        D = sp.diags(D, 0)
        wn = D @ w
    elif type == "gph":
        D = 1 / np.sqrt(D)
        D = sp.diags(D, 0)
        wn = D @ (w @ D)
    return wn


def Network_Enhancement(W_in: np.ndarray, order: int = 2, K: int = None, alpha: float = 0.9) -> np.ndarray:
    """
    Perform network enhancement.

    Args:
        W_in: The input matrix.
        order: The order of enhancement.
        K: The number of nearest neighbors to consider.
        alpha: The alpha parameter.

    Returns:
        The enhanced network matrix.
    """
    if K is None:
        K = min(20, int(np.ceil(len(W_in) / 10)))
    if alpha is None:
        alpha = 0.9

    W_in1 = W_in * (1 - np.eye(len(W_in)))
    zeroindex = np.where(np.sum(np.abs(W_in1), axis=1) > 0)[0]
    W0 = W_in[zeroindex, :][:, zeroindex]
    W = dn(W0, "ave")
    W = (W + W.T) / 2

    DD = np.sum(np.abs(W0), axis=1)

    if len(np.unique(W)) == 2:
        P = W
    else:
        P = dominateset(np.abs(W), min(K, len(W) - 1)) * np.sign(W)
    P = P + np.eye(len(P)) + np.diag(np.sum(np.abs(P), axis=1))
    P = TransitionFields(P)
    U, d, _ = np.linalg.svd(P)
    d = np.real(d - np.finfo(float).eps)
    d = (1 - alpha) * d / (1 - alpha * d**order)
    D = np.diag(np.real(d))
    W = U @ D @ U.T

    W = W * (1 - np.eye(len(W))) / (1 - np.diag(W))
    D = sp.diags(DD, 0)
    W = D @ W
    W[W < 0] = 0
    W = (W + W.T) / 2
    W_out = np.zeros_like(W_in)
    W_out[np.ix_(zeroindex, zeroindex)] = W

    return W_out
