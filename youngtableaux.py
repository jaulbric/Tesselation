#!/usr/bin/env python
import numpy as np
import warnings

def _kmin(p):
    N = len(p) + 1
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p, dtype=int)
    j = np.linspace(1, N-1, num=N-1, dtype=int)
    k = np.array([n * np.sum(j * p, dtype=int) / N - np.sum((n - j[:n-1]) * p[-1:N - n - 1:-1], dtype=int) for n in range(1, N)])
    return np.max(k)

def nonzero_p(p):
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p)
    return np.where(p != 0)[0]

def fj(p, i, l1, l2):
    N = len(p) + 1
    f = np.zeros(N+1)
    f[1:-1] = p

    f[l1+1] -= 1
    f[l1 + i - N + l2 + 2] += 1
    f[l2 + 1] -= 1

    return f[1:-1]

def l_pair(p, i):
    l1 = 0
    l2 = 0
    N = len(p) + 1
    L = nonzero_p(p)
    if np.all(L > N - i - 1):
        l1 = -1
        l2 = np.min(L)
    else:
        nmin = len(L) - np.argmax(L[::-1] <= N - 1 - i) - 1
        l1 = L[nmin]
        if nmin == len(L) - 1:
            # l1 = L[nmin]
            l2 = N-1
        else:
            # l1 = L[nmin]
            l2 = L[nmin+1]
    return l1, l2

def Nality(p):
    N = len(p) + 1
    j = np.linspace(1, N-1, N-1, dtype = int)
    return np.sum(j * p) % N

def kmin(p):
    try:
        assert type(p) == np.ndarray
    except AssertionError:
        p = np.asarray(p)

    _N = len(p) + 1
    _i = _N - Nality(p)
    _fj = np.array([])
    if _i == _N:
        warnings.warn("p has N-ality 0, meaning Q_i must be trivial.")
        _fj = p
    else:
        _l1, _l2 = l_pair(p, _i)
        _fj = fj(p, _i, _l1, _l2)
    return _kmin(_fj)
