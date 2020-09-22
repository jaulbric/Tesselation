#!/usr/bin/env python
import numpy as np

def _kmin(p):
    """Minimum number of copies of the adjoint such that tensor product X G_N ... G_N contains a trivial subspace. Used internally.

    Parameters:
        p : array_like
            Dynkin label of irreducible representation X.

    Returns:
        kmin : int
            Minimum number of copies of the adjoint representation such that the operator is gauge invariant.
    """
    N = len(p) + 1
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p, dtype=int)
    j = np.linspace(1, N - 1, num=N - 1, dtype=int)
    k = np.array([n * np.sum(j * p, dtype=int) / N - np.sum((n - j[:n - 1])
                                                            * p[-1:N - n - 1:-1], dtype=int) for n in range(1, N)], dtype=int)
    return np.max(k)


def nonzero_p(p):
    """Index of non zero Dynkin labels.

    Parameters:
        p : array_like
            Dynkin label.

    Returns:
        L : array_like
            Indicies i such that p[i] != 0
    """
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p, dtype=int)
    return np.where(p != 0)[0]

def _fj(p, i, l1, l2):
    """Optimum irreducible representation in the tensor product of X times Q_i. Used internally.

    Parameters:
        p  : array_like
            Dynkin label of X.
        i  : int
            Index of only nonzero Dynkin coefficient in Q_i.
        l1 : int
            largest l such that p[l] != 0 and l <= N - i.
        l2 : int
            Smallest l such that p[l] != 0 and l > l1.

    Returns:
        fj : array_like
            Resulting Dynkin label for optimal irreducible representation.
    """
    N = len(p) + 1
    f = np.pad(p, 1, mode="constant", constant_values=1)

    f[l1] -= 1
    f[l1 + i - N + l2] += 1
    f[l2] -= 1

    return f[1:-1]

def l_pair(p, i):
    """Pair of l that results in the optimal irreducible representation in tensor product of X times Q_i.

    Parameters:
        p : array_like
            Dynkin label of X.
        i  : int
            Index of only nonzero Dynkin coefficient in Q_i.

    Returns:
        l1, l2 : int, int
            Pair of l such that p[l1] !=0, p[l2] != 0, l1 is largest l such that l <= N - i, and is smallest l such that l > l1.
    """
    N = len(p) + 1
    tildep = np.pad(p, 1, mode="constant", constant_values=1)
    L = nonzero_p(tildep)
    l1 = np.max(L[L <= N - i])
    l2 = np.min(L[L > l1])
    return l1, l2


def Nality(p):
    """N-ality of irreducible representation X. The number of boxes in the Young tableau for X mod N.

    Parameters:
        p : array_like
            Dynkin label of X.

    Returns:
        t : int
            The number of boxes in the Young tableau for X mod N.
    """
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p)
    N = len(p) + 1
    j = np.linspace(1, N - 1, N - 1, dtype=int)
    return np.sum(j * p) % N


def kmin(p):
    """Minimum number of copies of the adjoint such that tensor product X Q_i G_N ... G_N  contains a trivial subspace.

    Parameters:
        p : array_like
            Dynkin label of irreducible representation X.

    Returns:
        kmin : int
            Minimum number of copies of the adjoint representation such that the operator is gauge invariant.
    """
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p)

    _N = len(p) + 1
    _i = _N - Nality(p)
    _fj_arr = np.array([])
    if _i == _N:
        print("Warning: p has N-ality 0, meaning Q_i must be trivial.")
        _fj_arr = p
    else:
        _l1, _l2 = l_pair(p, _i)
        _fj_arr = _fj(p, _i, _l1, _l2)
    return _kmin(_fj_arr)


def fj(p):
    """Optimum irreducible representation in the tensor product of X times Q_i.

    Parameters:
        p  : array_like

    Returns:
        fj : array_like
            Resulting Dynkin label for optimal irreducible representation.
    """
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p)

    _N = len(p) + 1
    _i = _N - Nality(p)
    _fj_arr = np.array([])
    if _i == _N:
        print("Warning: p has N-ality 0, meaning Q_i must be trivial.")
        _fj_arr = p
    else:
        _l1, _l2 = l_pair(p, _i)
        _fj_arr = _fj(p, _i, _l1, _l2)
    return _fj_arr
