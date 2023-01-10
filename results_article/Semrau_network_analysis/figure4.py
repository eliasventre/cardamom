# Network validation from ChIP-seq
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from sklearn.metrics import roc_curve, precision_recall_curve, auc

path = './inferred_networks/'

# List of algorithms
algos = ['CARDAMOM', 'HARISSA', 'GENIE3', 'SINCERITIES']

# List of colors
cmap = plt.get_cmap('tab10')
# color = {'CARDAMOM': cmap(4), 'HARISSA': cmap(3),
#     'GENIE3': cmap(0), 'SINCERITIES': cmap(1)}
color = {'CARDAMOM': 'purple', 'HARISSA': 'orangered',
    'GENIE3': 'dodgerblue', 'SINCERITIES': 'orange'}

# Load the data
NETWORK = np.load(path+'REAL.npy')
HARISSA = np.abs(np.load(path+'HARISSA.npy'))
GENIE3 = np.abs(np.load(path+'GENIE3.npy')).T
SINCERITIES = np.abs(np.load(path+'SINCERITIES.npy'))
CARDAMOM = np.abs(np.load(path+'CARDAMOM.npy'))

# Vectorize scores: directed version
vTd, vCd, vHd, vGd, vSd = [], [], [], [], []
G = 42
for i in range(G):
    for j in range(1, G):
        if i != j and NETWORK[i, j] >= 0:
            vCd.append(CARDAMOM[i, j])
            vTd.append(NETWORK[i, j])
            vHd.append(HARISSA[i, j])
            vGd.append(GENIE3[i, j])
            vSd.append(SINCERITIES[i, j])
vCd = np.array(vCd)
vTd = np.array(vTd)
vHd = np.array(vHd)
vGd = np.array(vGd)
vSd = np.array(vSd)

# Figure
fig = plt.figure(figsize=(8.115,2.4))
grid = gs.GridSpec(1, 2, wspace=0.1, width_ratios=[1,0.35])
panelAB = grid[0,0].subgridspec(1, 2, wspace=0.3)
panelC = grid[0,1]

# Axis settings
def configure(ax):
    ax.set_aspect('equal')
    ax.tick_params(axis='both', labelsize=7, pad=2, length=3)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)

# Panel settings
x0, y0 = -26, -3
opt = {'xy': (0,1), 'xycoords': 'axes fraction', 'fontsize': 10.5,
    'textcoords': 'offset points', 'annotation_clip': False}

# Curve settings
optc = {'lw': 1}

# Legend settings
optl = {'loc': 'lower right', 'fontsize': 6.5, 'frameon': False,
    'borderaxespad': 0.2, 'handlelength': 1.4, 'borderpad': 0.2,
    'labelspacing': 0.4, 'handletextpad': 0.5}

# A. ROC curves
ax = plt.subplot(panelAB[0])
ax.annotate('A', xytext=(x0,y0), fontweight='bold', **opt)
configure(ax)
# CARDAMOM
xCd, yCd, tCd = roc_curve(vTd, vCd)
a = int(100*auc(xCd,yCd))/100
l = f'CARDAMOM ({a})'
ax.plot(xCd, yCd, c=color['CARDAMOM'], label=l, **optc)
# HARISSA
xHd, yHd, tHd = roc_curve(vTd, vHd)
a = int(100*auc(xHd,yHd))/100
l = f'HARISSA ({a})'
ax.plot(xHd, yHd, c=color['HARISSA'], label=l, **optc)
# GENIE3
xGd, yGd, tGd = roc_curve(vTd, vGd)
a = int(100*auc(xGd,yGd))/100
l = f'GENIE3 ({a})'
ax.plot(xGd, yGd, c=color['GENIE3'], label=l, **optc)
# SINCERITIES
xSd, ySd, tSd = roc_curve(vTd, vSd)
a = int(100*auc(xSd,ySd))/100
l = f'SINCERITIES ({a})'
ax.plot(xSd, ySd, c= color['SINCERITIES'], label=l, **optc)
# Random estimator
l = f'Random ({0.5})'
ax.plot([0,1], [0,1], ls='--', c='lightgray', zorder=0, label=l, **optc)
# Legend
ax.legend(**optl)
ax.set_ylabel('True positive rate', fontsize=7)
ax.set_xlabel('False positive rate', fontsize=7)

# B. PR curves
ax = plt.subplot(panelAB[1])
ax.annotate('B', xytext=(x0,y0), fontweight='bold', **opt)
configure(ax)
# CARDAMOM
yCd, xCd, tCd = precision_recall_curve(vTd, vCd)
a = int(100*auc(xCd,yCd))/100
l = f'CARDAMOM ({a})'
ax.plot(xCd, yCd, c=color['CARDAMOM'], label=l, **optc)
# HARISSA
yHd, xHd, tHd = precision_recall_curve(vTd, vHd)
a = int(100*auc(xHd,yHd))/100
l = f'HARISSA ({a})'
ax.plot(xHd, yHd, c=color['HARISSA'], label=l, **optc)
# GENIE3
yGd, xGd, tGd = precision_recall_curve(vTd, vGd)
a = int(100*auc(xGd,yGd))/100
l = f'GENIE3 ({a})'
ax.plot(xGd, yGd, c=color['GENIE3'], label=l, **optc)
# SINCERITIES
ySd, xSd, tSd = precision_recall_curve(vTd, vSd)
a = int(100*auc(xSd,ySd))/100
l = f'SINCERITIES ({a})'
ax.plot(xSd, ySd, c= color['SINCERITIES'], label=l, **optc)
# Random estimator
m = np.mean(vTd)
a = int(100*m)/100
l = f'Random ({a})'
ax.plot([0,1], [m,m], ls='--', c='lightgray', zorder=0, label=l, **optc)
# Legend
ax.legend(**optl)
ax.set_ylabel('Precision', fontsize=7)
ax.set_xlabel('Recall', fontsize=7)


# Export the figure
fig.savefig('figure4.pdf', bbox_inches='tight', pad_inches=0.04)
