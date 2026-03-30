import numpy as np
from scipy import sparse

def _get_mean_var(X):
    if sparse.issparse(X):
        mean = np.array(X.mean(axis=0)).ravel()
        var = np.array(X.power(2).mean(axis=0)).ravel() - mean**2
    else:
        mean = np.mean(X, axis=0)
        var = np.var(X, axis=0, ddof=1)
    return mean, var

