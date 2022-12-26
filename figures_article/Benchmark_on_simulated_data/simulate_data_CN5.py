# Generate data for the 4-gene network of figure 3
import sys; sys.path += ['../']
import numpy as np
from harissa import NetworkModel

np.random.seed(0)

# Number of cells
C = 1000

# Number of runs
N = 10

# Time points
t = [0, 6, 12, 24, 36, 48, 60, 72, 84, 96] # np.linspace(0, 100, 10)
k = np.linspace(0, C, len(t) + 1, dtype='int')
print(f't = {t}')
time = np.zeros(C, dtype='int')
for i in range(len(t)): time[k[i]:k[i+1]] = t[i]

# Number of genes
G = 5

# Prepare data
data = np.zeros((C+1,G+2), dtype='int')
data[0][1:] = np.arange(G+1)
data[1:,0] = time # Time points
data[1:,1] = 100 * (time > 0) # Stimulus

for r in range(N):
    print(f'Run {r+1}...')
    
    # Initialize the model
    model = NetworkModel(G)
    model.d[0] = 0.5
    model.d[1] = 0.1

    model.basal[1:] = [-5, 4, 4, -5, -5]
    model.inter[0, 1] = 10
    model.inter[1, 2] = -10
    model.inter[2, 3] = -10
    model.inter[3, 4] = 10
    model.inter[4, 5] = 10
    model.inter[5, 1] = -10

    # Save true network topology
    inter = 1 * (abs(model.inter) > 0)
    np.save(f'CN5/True/inter_{r+1}', inter)

    # Generate data
    for k in range(C):
        # print(f'* Cell {k+1}')
        sim = model.simulate(time[k], burnin=5)
        data[k+1,2:] = np.random.poisson(sim.m[-1])

    # Save data for use with PIDC
    fname = f'CN5/Data/data_{r+1}.txt'
    np.savetxt(fname, data.T, fmt='%d', delimiter='\t')
