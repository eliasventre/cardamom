# Script pour lancer tous les benchmarks de HARISSA
import sys; sys.path += ['../']
import time as timer
import numpy as np
from harissa import NetworkModel

# Number of runs
N = 10

# Inference for Cycle
for r in range(N):
    fname = 'CN5/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
    x = data.T
    x[:,0] = time
    model = NetworkModel()
    model.fit(x)
    score = model.inter
    np.save('CN5/HARISSA/score_{}'.format(r+1), score)

# Inference for Network4
for r in range(N):
    fname = 'FN4/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
    x = data.T
    x[:,0] = time
    model = NetworkModel()
    model.fit(x)
    score = model.inter
    np.save('FN4/HARISSA/score_{}'.format(r+1), score)


# Inference for Trifurcation
for r in range(N):
    fname = 'FN8/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
    x = data.T
    x[:,0] = time
    model = NetworkModel()
    model.fit(x)
    score = model.inter
    np.save('FN8/HARISSA/score_{}'.format(r+1), score)

# Inference for Bifurcation
for r in range(N):
    fname = 'BN8/Data/data_{}.txt'.format(r+1)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
    time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
    x = data.T
    x[:,0] = time
    model = NetworkModel()
    model.fit(x)
    score = model.inter
    np.save('BN8/HARISSA/score_{}'.format(r+1), score)

# Inference for tree-like networks
for n in [5,10,20,50,100]:
    runtime = np.zeros(N)
    for r in range(N):
        fname = 'Trees{}/Data/data_{}.txt'.format(n,r+1)
        data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:,1:]
        time = np.loadtxt(fname, dtype=int, delimiter='\t')[0,1:]
        x = data.T
        x[:,0] = time
        model = NetworkModel()
        t0 = timer.time()
        model.fit(x)
        t1 = timer.time()
        runtime[r] = t1 - t0
        score = model.inter
        np.save('Trees{}/HARISSA/score_{}'.format(n,r+1), score)
    # Save running times
    np.savetxt('Trees{}/HARISSA/runtime.txt'.format(n), runtime.T)


