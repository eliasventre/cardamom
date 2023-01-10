# Benchmark (directed|undirected) for all test networks
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from sklearn.metrics import precision_recall_curve, auc

# Data location
path = './'

# Directed inference from snapshots
algoD = ['HARISSA', 'CARDAMOM', 'GENIE3', 'SINCERITIES']

# Undirected inference from snapshots
algoU = ['HARISSA', 'CARDAMOM', 'GENIE3', 'SINCERITIES', 'PIDC', 'PEARSON']


# List of benchmarks
benchmarks = ['FN4', 'CN5', 'FN8', 'BN8',
    'Trees5', 'Trees10', 'Trees20', 'Trees50', 'Trees100']

# Number of datasets for each benchmark
N = 10

# List of colors
cmap = plt.get_cmap('tab20')
c = {'HARISSA': (cmap(6), cmap(7)), 'CARDAMOM': (cmap(8), cmap(9)),
    'GENIE3': (cmap(0), cmap(1)), 'SINCERITIES': (cmap(2), cmap(3)),
    'PIDC': (cmap(4), cmap(5)), 'PEARSON': (cmap(14), cmap(15)),
    'SCRIBE_timed': (cmap(10), cmap(11)), 'SCRIBE_wadd': (cmap(12), cmap(13)),
    'SCRIBE_pseudotimed': (cmap(16), cmap(17)), 'Random': 2*('lightgray',)}

# Figure
fig = plt.figure(figsize=(7.85,8))
grid = gs.GridSpec(5, 3, hspace=0.3, wspace=0.18,
    height_ratios=[1,1,1,1,1], width_ratios=[3.7,2.5,5.6])

# Axis settings
def configure(ax):
    w = 0.7
    ax.tick_params(direction='out', length=3, width=w)
    ax.tick_params(axis='x', pad=2, labelsize=5.5)
    ax.tick_params(axis='y', pad=0.5, labelsize=5.5)
    for x in ['top','bottom','left','right']: ax.spines[x].set_linewidth(w)
    ax.set_ylim(0,1)

# Boxplot settings
opt_box = {'patch_artist': True, 'widths': [.25]}
def configure_box(box, c):
    w = 0.8
    for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(box[item], color=c[0], lw=w)
    plt.setp(box['boxes'], facecolor=c[1])
    plt.setp(box['fliers'], markeredgecolor=c[0], ms=3, markerfacecolor=c[1],
        markeredgewidth=w)

# Panel settings
x, y = -11, 10
xn, yn = -0.142, 0.875
opt = {'xy': (0,1), 'xycoords': 'axes fraction', 'fontsize': 10,
    'textcoords': 'offset points', 'annotation_clip': False}

# Initialize data for trees
auprTreesD = {algo: [] for algo in algoD+['Random']}
auprTreesU = {algo: [] for algo in algoU+['Random']}

# Load the data
for n, benchmark in enumerate(benchmarks):
    # Initialize values
    auprD = {algo: [] for algo in algoD+['Random']}
    auprU = {algo: [] for algo in algoU+['Random']}

    # Routine for each dataset
    for r in range(1,N+1):
        # True network
        inter = abs(np.load(path+f'{benchmark}/True/inter_{r}.npy'))
        G = inter.shape[0]

        # 1. Directed inference from snapshots
        edges = [(i,j) for i in range(G) for j in set(range(1,G))-{i}]
        y0 = np.array([inter[i,j] for (i,j) in edges])
        auprD['Random'].append(np.mean(y0))
        for algo in algoD:
            score = abs(np.load(path+f'{benchmark}/{algo}/score_{r}.npy'))
            if algo=='GENIE3': score = score.T # Solve direction problem?
            y1 = np.array([score[i,j] for (i,j) in edges])
            precision, recall, thresholds = precision_recall_curve(y0, y1)
            auprD[algo].append(auc(recall,precision))

        # 2. Undirected inference from snapshots
        edges = [(i,j) for i in range(G) for j in range(i+1,G)]
        y0 = np.array([max(inter[i,j],inter[j,i]) for (i,j) in edges])
        auprU['Random'].append(np.mean(y0))
        for algo in algoU:
            score = abs(np.load(path+f'{benchmark}/{algo}/score_{r}.npy'))
            y1 = np.array([max(score[i,j],score[j,i]) for (i,j) in edges])
            precision, recall, thresholds = precision_recall_curve(y0, y1)
            auprU[algo].append(auc(recall,precision))


    # Store data for trees
    if benchmark[:5] == 'Trees':
        for algo in algoD+['Random']:
            auprTreesD[algo].append(np.mean(auprD[algo]))
        for algo in algoU+['Random']:
            auprTreesU[algo].append(np.mean(auprU[algo]))

    # Non-tree benchmarks
    name = ['FN4', 'CN5', 'FN8', 'BN8']
    if n < 4:
        # A. Snapshot-based
        ax = plt.subplot(grid[n,0])
        configure(ax)
        if n == 0:
            title = 'Snapshot-based'
            ax.annotate('A', xytext=(x+0.2,y), fontweight='bold', **opt)
            ax.annotate(title, xytext=(x+14,y), **opt)
        # Draw the baseline
        b = np.mean(auprD['Random'])
        ax.plot([0,5], [b,b], color='lightgray', ls='--', lw=0.8, zorder=0)
        ax.set_xlim(0.8,4.5)
        # Draw the histograms
        for i, algo in enumerate(algoD):
            box = ax.boxplot([auprD[algo]], positions=[i+1], **opt_box)
            configure_box(box, c[algo])
        ax.set_xticklabels(algoD)
        ax.set_ylabel('AUPR', fontsize=6)
        # Network name
        optn = {'fontsize': 9, 'transform': ax.transAxes, 'ha': 'right'}
        ax.text(xn, yn, name[n], **optn)
        ax.text(xn, yn+0.01, name[n], color='none', zorder=0, bbox=dict(
            boxstyle='round,pad=0.2',fc='none',ec='lightgray',lw=0.8), **optn)


        # C. Snapshot-based (undirected edges)
        ax = plt.subplot(grid[n,2])
        configure(ax)
        if n == 0:
            title = 'Snapshot-based (undirected edges)'
            ax.annotate('B', xytext=(x,y), fontweight='bold', **opt)
            ax.annotate(title, xytext=(x+14,y), **opt)
        # Draw the baseline
        b = np.mean(auprU['Random'])
        ax.plot([0,7], [b,b], color='lightgray', ls='--', lw=0.8, zorder=0)
        ax.set_xlim(0.8,6.4)
        # Draw the histograms
        for i, algo in enumerate(algoU):
            box = ax.boxplot([auprU[algo]], positions=[i+1], **opt_box)
            configure_box(box, c[algo])
        ax.set_xticklabels(algoU)
        # ax.tick_params(axis='y', labelleft=False)

# Tree benchmarks
s = {'ls': '--', 'lw': 0.85, 'marker': '.', 'ms': 4}
p = {'borderaxespad': 0, 'frameon': False, 'fontsize': 5.5,
    'handlelength': 1.2, 'handletextpad': 0.5}
size = [5, 10, 20, 50, 100]

# A. Snapshot-based
ax = plt.subplot(grid[4,0])
configure(ax)
# Draw the curves
for algo in algoD+['Random']:
    ax.plot(size, auprTreesD[algo], color=c[algo][0], label=algo, **s)
ax.legend(loc='upper right', **p)
ax.set_xticks(size)
ax.set_xlabel('No. of genes', fontsize=6, labelpad=1.8)
ax.set_ylabel('AUPR', fontsize=6)
# Network name
optn = {'fontsize': 9, 'transform': ax.transAxes, 'ha': 'right'}
ax.text(xn, yn, 'Trees', **optn)
ax.text(xn, yn+0.01, 'Trees', color='none', zorder=0, bbox=dict(
    boxstyle='round,pad=0.2',fc='none',ec='lightgray',lw=0.8), **optn)


# C. Snapshot-based (undirected)
ax = plt.subplot(grid[4,2])
configure(ax)
# Draw the curves
for algo in algoU+['Random']:
    ax.plot(size, auprTreesU[algo], color=c[algo][0], label=algo, **s)
ax.legend(loc='upper right', ncol=2, **p)
ax.set_xticks(size)
ax.set_xlabel('No. of genes', fontsize=6, labelpad=1.8)

# Export the figure
fig.savefig('figure2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.04)
