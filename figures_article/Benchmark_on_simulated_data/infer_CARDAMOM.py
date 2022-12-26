# Script pour lancer tous les benchmarks de CARDAMOM
import sys; sys.path += ['./../../']
import time as timer
import numpy as np
from cardamom import NetworkModel

# Number of runs
N = 10
# Print information
verb = 0

# Inference for Network4
for r in range(0, N):
    fname = 'FN4/Data/data_{}.txt'.format(r + 1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
    x = data.T
    x[:, 0] = time
    G = np.size(x, 1)
    model = NetworkModel(G - 1)
    model.fit(x, verb=verb)
    np.save('FN4/CARDAMOM/score_{}'.format(r+1), model.inter)

# Inference for Cycle
for r in range(0, N):
    fname = 'CN5/Data/data_{}.txt'.format(r + 1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:, 1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0, 1:]
    x = data.T
    x[:, 0] = time
    G = np.size(x, 1)
    model = NetworkModel(G - 1)
    model.fit(x, verb=verb)
    np.save('CN5/CARDAMOM/score_{}'.format(r + 1), model.inter)

# Inference for Bifurcation
for r in range(0, N):
    fname = 'BN8/Data/data_{}.txt'.format(r + 1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:, 1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0, 1:]
    x = data.T
    x[:, 0] = time
    G = np.size(x, 1)
    model = NetworkModel(G - 1)
    model.fit(x, verb=verb)
    np.save('BN8/CARDAMOM/score_{}'.format(r + 1), model.inter)


# Inference for Trifurcation
for r in range(0, N):
    fname = 'FN8/Data/data_{}.txt'.format(r + 1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:, 1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0, 1:]
    x = data.T
    x[:, 0] = time
    G = np.size(x, 1)
    model = NetworkModel(G - 1)
    model.fit(x, verb=verb)
    np.save('FN8/CARDAMOM/score_{}'.format(r + 1), model.inter)

# Inference for tree-like networks
for n in [5, 10, 20, 50, 100]:
    runtime = np.zeros(N)
    for r in range(N):
        fname = 'Trees{}/Data/data_{}.txt'.format(n,r+1)
        data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
        time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
        x = data.T
        x[:,0] = time
        G = np.size(x, 1)
        model = NetworkModel(G - 1)
        t0 = timer.time()
        model.fit(x, verb=verb)
        t1 = timer.time()
        runtime[r] = t1 - t0
        score = model.inter
        np.save('Trees{}/CARDAMOM/score_{}'.format(n,r+1), score)
    # Save running times
    np.savetxt('Trees{}/CARDAMOM/runtime.txt'.format(n), runtime.T)
