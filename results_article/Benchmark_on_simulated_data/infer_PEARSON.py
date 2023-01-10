# Script pour lancer tous les benchmarks de PEARSON
import sys; sys.path += ['./_scripts']
import time as timer
import numpy as np
from scipy import stats

# Number of runs
N = 10

# Inference for Network4
for r in range(N):
    fname = 'FN4/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    x = data.T
    G = np.size(x, 1)
    score = np.zeros((G, G))
    for i in range(0, G):
        for j in range(0, G):
            score[i, j] = stats.pearsonr(x[:, i], x[:, j])[0]
    np.save('FN4/PEARSON/score_{}'.format(r+1), score)

# Inference for Cycle
for r in range(N):
    fname = 'CN5/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    x = data.T
    G = np.size(x, 1)
    score = np.zeros((G, G))
    for i in range(0, G):
        for j in range(0, G):
            score[i, j] = stats.pearsonr(x[:, i], x[:, j])[0]
    np.save('CN5/PEARSON/score_{}'.format(r+1), score)

# Inference for Bifurcation
for r in range(N):
    fname = 'BN8/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    x = data.T
    G = np.size(x, 1)
    score = np.zeros((G, G))
    for i in range(0, G):
        for j in range(0, G):
            score[i, j] = stats.pearsonr(x[:, i], x[:, j])[0]
    np.save('BN8/PEARSON/score_{}'.format(r+1), score)


# Inference for Trifurcation
for r in range(N):
    fname = 'FN8/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    x = data.T
    G = np.size(x, 1)
    score = np.zeros((G, G))
    for i in range(0, G):
        for j in range(0, G):
            score[i, j] = stats.pearsonr(x[:, i], x[:, j])[0]
    np.save('FN8/PEARSON/score_{}'.format(r+1), score)

# Inference for tree-like networks
for n in [5, 10, 20,50,100]:
    runtime = np.zeros(N)
    for r in range(N):
        fname = 'Trees{}/Data/data_{}.txt'.format(n,r+1)
        data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
        x = data.T
        t0 = timer.time()
        G = np.size(x, 1)
        score = np.zeros((G, G))
        for i in range(0, G):
            for j in range(0, G):
                score[i, j] = stats.pearsonr(x[:, i], x[:, j])[0]
        t1 = timer.time()
        runtime[r] = t1 - t0
        np.save('Trees{}/PEARSON/score_{}'.format(n,r+1), score)
    # Save running times
    np.savetxt('Trees{}/PEARSON/runtime.txt'.format(n), runtime.T)

