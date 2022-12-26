# Network inferred by CARDAMOM
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import sys; sys.path += ['../..']
from matplotlib.lines import Line2D
from utils import plot_network as nt
nt.inhib = '#F03E3E'

# Load the data
path = './'
names = np.loadtxt(path+'../../Semrau/Data/panel_genes.txt', dtype='str')[:,1]
time = np.array([0, 6, 12, 24, 36, 48, 60, 72, 96])
theta_real = np.load(path+'inferred_networks/REAL.npy')
inter = np.load(path+'inferred_networks/CARDAMOM.npy') * 10

# Node positions
p = np.load('utils/positions_edges.npy')

# Global transform
px = 1*p
px[:,0] += 0.75
px[:,1] += 0.25
scale = 1.32

# Network links
G = inter.shape[0]
theta = inter * (1-np.eye(G))
cnt = 0
seuil_min, seuil_max = 4, 4.5
theta[0,:] *= (np.abs(theta[0,:]) >= 4)
for j in range(1,G):
    for i in range(1,G):
        # Negative threshold
        if len(theta[theta[:,i]<0,i]):
            seuil_neg = np.quantile(-theta[theta[:,i]<0,i], 0.96)
        else: seuil_neg = seuil_min
        # Positive threshold
        if len(theta[theta[:,i]>0,i]):
            seuil_pos = np.quantile(theta[theta[:,i]>0,i], 0.96)
        else: seuil_pos = seuil_min
        # Selection of edges
        s0 = min(max(seuil_neg,seuil_min),seuil_max)
        s1 = min(max(seuil_pos,seuil_min),seuil_max)
        if theta[j,i] < -s0 or theta[j,i] > s1: cnt += 1
        else: theta[j,i] = 0

# Gene group colors
# c1, c2, c3, c4 = ['blue', 'dodgerblue', 'brown', 'orange']
c1, c2, c4 = [plt.get_cmap('tab10')(i) for i in [0,4,1]]
c3 = plt.get_cmap('Paired')(11)
vcolor = ['gray'] + 10*[c1] + 10*[c2] + 11*[c3] + 10*[c4]

# Figure
fig = plt.figure(figsize=(8.023,8.075))
grid = gs.GridSpec(1, 1)
ax = plt.subplot(grid[0,0])

# Plot the network
nt.plot_network(theta, axes=ax, scale=scale, layout=px, names=names,
    vcolor=vcolor, fontsize=9, nodesize=1.4, alpha=0.5)

# Add known interactions
pos = ax.get_position()
scale *= np.min([pos.width,pos.height])
for k1, k2 in zip(*theta.nonzero()):
    if theta_real[k1,k2] == 1: weight_real = 1
    elif theta_real[k1,k2] == 0: weight_real = -1
    else: weight_real = 0
    if weight_real != 0:
        fill = False
        if weight_real > 0: fill = True
        s = 0.37
        if (k1,k2) == (0,4): s = 0.35
        if (k1,k2) == (0,13): s = 0.35
        if (k1,k2) == (0,33): s = 0.45
        if (k1,k2) == (1,2): s = 0.45
        if (k1,k2) == (1,3): s = 0.40
        if (k1,k2) == (1,34): s = 0.515
        if (k1,k2) == (37,32): s = 0.2
        if (k1,k2) == (2,13): s = 0.17
        if (k1,k2) == (2,31): s = 0.11
        if (k1,k2) == (6,15): s = 0.34
        x1, y1 = scale*px[k1]
        x2, y2 = scale*px[k2]
        circle = plt.Circle((x1 + s*(x2-x1), y1 + s*(y2-y1)),
            radius=0.024, fill=fill, color='k', lw=0.8)
        ax.add_artist(circle)

ax.axis('off')
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
# print(np.sum(theta_real[0]==1))

# Gene groups
labels = ['Pluripotency', 'Post-implantation epiblast     ',
    'Extraembryonic endoderm', 'Neuroectoderm']
lines = [Line2D([0],[0],color=c,lw=4) for c in [c1,c2,c4,c3]]
ax.legend(lines, labels, ncol=2, frameon=False, columnspacing=2,
    labelspacing=0.5, bbox_to_anchor=(0.15,1), handletextpad=0.8,
    loc='upper left', handlelength=1, fontsize=11)

# Export the figure
fig.savefig('figure5.pdf', bbox_inches='tight', pad_inches=0.04)
