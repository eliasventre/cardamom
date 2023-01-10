"""
Inference of basal parameters
This file comes from the 'Harissa' package written by Ulysse Herbach
        https://github.com/ulysseherbach/harissa
"""
import numpy as np
from numpy import log
from scipy.special import psi, polygamma, gammaln

def log_gamma_poisson_pdf(k, a, b):

    ga, gk = gammaln(a), gammaln(k + 1)
    return gammaln(a + k) - (ga + gk) + a * np.log(b) - (a + k) * np.log(b + 1)

def estim_gamma_poisson(x):
    """
    Estimate parameters a and b of the Gamma-Poisson(a,b) distribution,
    a.k.a. negative binomial distribution, using the method of moments.
    """
    m1 = np.mean(x)
    m2 = np.mean(x*(x-1))
    if m1 == 0: return 0, 1
    r = m2 - m1**2
    if r > 0: a = m1 ** 2 / r
    else:
        v = np.var(x)
        if v == 0: return 0, 1
        a = m1 ** 2 /v
    b = a / m1
    return a, b

def infer_kinetics(x, times, tol=1e-5, max_iter=1000, verb=False):
    """
    Infer parameters a[0], ..., a[m-1] and b of a Gamma-Poisson model
    with time-dependant a and constant b for a given gene at m time points.

    Parameters
    ----------
    x[k] = gene expression in cell k
    times[k] = time point of cell k
    """

    t = np.sort(list(set(times)))
    m = t.size
    n = np.zeros(m) # Number of cells for each time point
    a = np.zeros(m)
    b = np.zeros(m)
    # Initialization of a and b
    for i in range(m):
        cells = (times == t[i])
        n[i] = np.sum(cells)
        a[i], b[i] = estim_gamma_poisson(x[cells])
    b = np.mean(b)
    # Newton-like method
    k, c = 0, 0
    sx = np.sum(x)
    while (k == 0) or (k < max_iter and c > tol):
        da = np.zeros(m)
        for i in range(m):
            if a[i] > 0:
                cells = (times == t[i])
                z = a[i] + x[cells]
                p0 = np.sum(psi(z))
                p1 = np.sum(polygamma(1, z))
                d = n[i]*(log(b)-log(b+1)-psi(a[i])) + p0
                h = p1 - n[i]*polygamma(1, a[i])
                da[i] = -d/h
        anew = a + da
        if np.sum(anew < 0) == 0: a[:] = anew
        else:
            max_test = 5
            test = 0
            da *= 0.5
            while (np.sum(a + da < 0) > 0) and (test < max_test):
                da *= 0.5
                test += 1
            if test < max_test: a[:] = a + da
            else: print('Warning: parameter a not improved')
        if np.sum(a == 0) == 0:
            b = np.sum(n*a)/sx
        else: b = b
        c = np.max(np.abs(da))
        k += 1
    if (k == max_iter) and (c > tol):
        print('Warning: bad convergence (b = {})'.format(b))
        a, b = a/b, 1
    # if verb: print('Estimation done in {} iterations'.format(k))
    if np.sum(a < 0) > 0: print('WARNING: a < 0')
    if b < 0: print('WARNING: b < 0')
    if np.all(a == 0): print('WARNING: a == 0')
    # if k > 20 and np.max(a/b) > 2: print(k, np.max(a/b))

    return a, b
