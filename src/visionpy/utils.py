import numpy as np
from scipy import sparse

def _get_mean_var(X, axis=0):
    if sparse.issparse(X):
        mean = np.array(X.mean(axis=axis)).ravel()
        var = np.array(X.power(2).mean(axis=axis)).ravel() - mean**2
    else:
        mean = np.mean(X, axis=axis)
        var = np.var(X, axis=axis, ddof=1)
    return mean, var

