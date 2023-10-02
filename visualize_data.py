### PLOT THE UMAP

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import ks_2samp as ks
from umap import UMAP
import sys, getopt


plot_distributions = 1 # for plotting the marginal distributions of the simulated genes
plot_umap = 1 # for plotting the UMAP reduction of the simulated dataset
plot_pvalues = 0 # for plotting the comparison of the marginals using a Kolmogorov-Smornov test
verb = 1


def configure(ax):
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('UMAP1', fontsize=7, weight='bold')
    ax.set_ylabel('UMAP2', fontsize=7, weight='bold')
    ax.tick_params(axis='both', labelsize=7)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines[['top', 'right']].set_visible(False)


def plot_data_distrib(data_reference, data_simulated, t_real, t_netw, names, file):

    if len(t_real) != len(t_netw):
        print('Times are not the same !')
        return None
    rat = 5
    nb_by_pages = 10
    nb_genes = len(names)
    list_genes = np.arange(nb_genes)
    nb_pages = int(nb_genes / nb_by_pages) + 1
    with PdfPages('./{}/Results/Marginals.pdf'.format(file)) as pdf:
        for i in range(nb_pages):
            fig, ax = plt.subplots(len(t_netw), min(nb_by_pages, nb_genes),
                                   figsize=(min(nb_by_pages, nb_genes) * rat, len(t_netw) * rat))
            if nb_genes - i*nb_by_pages < nb_by_pages and nb_by_pages < nb_genes:
                for j in range(nb_genes - i*nb_by_pages, nb_by_pages):
                    for cnt_t, time in enumerate(t_real):
                        ax[cnt_t, j].set_axis_off()
            for cnt_g, g in enumerate(list_genes[i*nb_by_pages:min((i+1)*nb_by_pages, nb_genes)]):
                n_max = max(np.quantile(data_reference[g, :], 1), np.quantile(data_simulated[g, :], 1)) + 1
                n_bins = min(int(n_max / 2) + 1, 25)
                for cnt_t, time in enumerate(t_real):
                    data_tmp_simulated = data_simulated[g, data_simulated[0, :] == t_netw[cnt_t]]
                    data_tmp_reference = data_reference[g, data_reference[0, :] == t_real[cnt_t]]
                    if time == t_netw[-1]: ax[-1, cnt_g].set_xlabel('mRNA (copies per cell)', fontsize=20)
                    if time == t_netw[0]: ax[cnt_t, cnt_g].set_title(names[g], fontweight="bold", fontsize=30)
                    ax[cnt_t, cnt_g].hist(data_tmp_reference, density=True, bins=np.linspace(0, n_max, n_bins),
                                        color='grey', histtype='bar', alpha=0.7)
                    ax[cnt_t, cnt_g].hist(data_tmp_simulated, density=True, bins=np.linspace(0, n_max, n_bins),
                                                  ec='red', histtype=u'step', alpha=1, linewidth=2)
                    ax[cnt_t, cnt_g].legend(labels=['Model (t = {}h)'.format(int(t_real[cnt_t])),
                                                  'Data (t = {}h)'.format(int(t_real[cnt_t]))])
            pdf.savefig(fig)
            plt.close()

def plot_data_umap(data_real, data_netw, t_real, t_netw, inputfile):
    data_real = data_real[:, :]
    data_netw = data_netw[:, :]

    # Compute the UMAP projection
    reducer = UMAP(random_state=42, min_dist=0.7)
    proj = reducer.fit(data_real[1:,:].T)
    x_real = proj.transform(data_real[1:,:].T)
    x_netw = proj.transform(data_netw[1:,:].T)

    # Figure
    fig = plt.figure(figsize=(10, 4))
    grid = gs.GridSpec(2, 2, height_ratios=[1, 0.05], wspace=0.3)
    ax0 = plt.subplot(grid[0, 0])
    ax1 = plt.subplot(grid[0, 1])
    ax3 = plt.subplot(grid[1, :])

    # Panel settings
    opt = {'xy': (0, 1), 'xycoords': 'axes fraction', 'fontsize': 10,
           'textcoords': 'offset points', 'annotation_clip': False}

    # Timepoint colors
    T = len(t_real)
    cmap = [plt.get_cmap('viridis', T)(i) for i in range(T)]
    c_real = [cmap[np.argwhere(t_real==t)[0,0]] for t in data_real[0, :]]
    c_netw = [cmap[np.argwhere(t_netw==t)[0,0]] for t in data_netw[0, :]]

    # A. Original data
    configure(ax0)
    title = 'Original data'
    ax0.annotate('A', xytext=(-11, 6), fontweight='bold', **opt)
    ax0.annotate(title, xytext=(3, 6), **opt)
    ax0.scatter(x_real[:, 0], x_real[:, 1], c=c_real, s=2)

    # B. Inferred network
    configure(ax1)
    title = 'Inferred network'
    ax1.annotate('B', xytext=(-11, 6), fontweight='bold', **opt)
    ax1.annotate(title, xytext=(3, 6), **opt)
    ax1.scatter(x_netw[:, 0], x_netw[:, 1], c=c_netw, s=2)
    ax1.set(xlim=ax0.get_xlim(), ylim=ax0.get_ylim())

    # Legend panel
    labels = [f'{int(t_real[k])}h' for k in range(T)]
    lines = [Line2D([0], [0], color=cmap[k], lw=5) for k in range(T)]
    ax3.legend(lines, labels, ncol=T, frameon=False, borderaxespad=0,
               loc='lower right', handlelength=1, fontsize=8.5)
    ax3.text(-0.02, 0.8, 'Timepoints:', transform=ax3.transAxes, fontsize=8.5)
    ax3.axis('off')

    # Export the figure
    fig.savefig('./{}/Results/UMAP.pdf'.format(inputfile), dpi=300, bbox_inches='tight', pad_inches=0.02)

def compare_marginals(data_real, data_netw, t_real, t_netw, genes, file):


    T = len(t_real)
    G = len(genes)

    pval_netw = np.ones((T, G))
    for cnt_t in range(T):
        data_tmp_real = data_real[:,data_real[0,:] == t_real[cnt_t]]
        data_tmp_netw = data_netw[:,data_netw[0,:] == t_netw[cnt_t]]
        for cnt_g in range(0,G):
            stat_tmp = ks(data_tmp_real[cnt_g, :], data_tmp_netw[cnt_g, :])
            pval_netw[cnt_t, cnt_g] = stat_tmp[1]

    # Figure
    fig = plt.figure(figsize=(8*G/42,8.1*G/42))
    grid = gs.GridSpec(6, 4, wspace=0, hspace=0,
        width_ratios=[0.09,1.48,0.32,1],
        height_ratios=[0.49,0.2,0.031,0.85,0.22,0.516])
    panelA = grid[0,:]
    # Panel settings
    opt = {'xy': (0,1), 'xycoords': 'axes fraction', 'fontsize': 10,
        'textcoords': 'offset points', 'annotation_clip': False}

    # Color settings
    colors = ['#d73027','#f46d43','#fdae61','#fee08b','#ffffbf',
        '#d9ef8b','#a6d96a','#66bd63','#1a9850']

    # A. KS test p-values
    axA = plt.subplot(panelA)
    axA.annotate('A', xytext=(-14,6), fontweight='bold', **opt)
    axA.annotate('KS test p-values', xytext=(0,6), **opt)
    # axA.set_title('KS test p-values', fontsize=10)
    cmap = LinearSegmentedColormap.from_list('pvalue', colors)
    norm = Normalize(vmin=0, vmax=0.1)
    # Plot the heatmap
    im = axA.imshow(pval_netw, cmap=cmap, norm=norm)
    axA.set_aspect('equal','box')
    axA.set_xlim(-0.5,G-0.5)
    axA.set_ylim(T-0.5,-0.5)
    # Create colorbar
    divider = make_axes_locatable(axA)
    cax = divider.append_axes('right', '1.5%', pad='2%')
    cbar = axA.figure.colorbar(im, cax=cax, extend='max')
    pticks = np.array([0,1,3,5,7,9])
    cbar.set_ticks(pticks/100 + 0.0007)
    cbar.ax.set_yticklabels([0]+[f'{p}%' for p in pticks[1:]], fontsize=6)
    cbar.ax.spines[:].set_visible(False)
    cbar.ax.tick_params(axis='y',direction='out', length=1.5, pad=1.5)
    axA.set_xticks(np.arange(G))
    axA.set_yticks(np.arange(T))
    axA.set_xticklabels(genes, rotation=45, ha='right', rotation_mode='anchor',
        fontsize=3)
    axA.set_yticklabels([f'{int(t)}h' for t in t_real], fontsize=6.5)
    axA.spines[:].set_visible(False)
    axA.set_xticks(np.arange(G+1)-0.5, minor=True)
    axA.set_yticks(np.arange(T+1)-0.5, minor=True)
    axA.grid(which='minor', color='w', linestyle='-', linewidth=1)
    axA.tick_params(which='minor', bottom=False, left=False)
    axA.tick_params(which='major', bottom=False, left=False)
    axA.tick_params(axis='x',direction='out', pad=-0.1)
    axA.tick_params(axis='y',direction='out', pad=-0.1)

    # Export the figure
    fig.savefig('./{}/Results/Comparison.pdf'.format(file), dpi=300, bbox_inches='tight', pad_inches=0.02)

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:", ["ifile="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg

    ### PLOT DISTRIBUTION

    p = '{}/'.format(inputfile)  # Name of the file where are the data

    # Load  data
    # Data location
    path_real = p + 'Data/panel_real.txt'
    path_netw = p + 'Data/panel_simulated.txt'

    # Load the data
    data_real = np.loadtxt(path_real, dtype=float, delimiter='\t')[1:, 1:]
    data_real[0, :] = np.loadtxt(path_real, dtype=float, delimiter='\t')[0, 1:]
    data_netw = np.loadtxt(path_netw, dtype=float, delimiter='\t')[1:, 1:]
    data_netw[0, :] = np.loadtxt(path_netw, dtype=float, delimiter='\t')[0, 1:]

    # Names of genes
    names = np.loadtxt(p + 'Data/panel_genes.txt', dtype='str')[:, 1]  # Load  the names of the genes
    t_real = list(set(data_real[0, :]))  # Time of each sample for the reference datasets
    t_real.sort()
    t_netw = list(set(data_netw[0, :]))
    t_netw.sort()

    # Sort genes
    order = np.argsort(names[1:])
    names = [names[0]] + list(names[order+1])

    # Order the data
    order = [0] + list(order+1)
    data_real = data_real[order,:]
    data_netw = data_netw[order,:]

    if plot_distributions:
        plot_data_distrib(data_real, data_netw, t_real, t_netw, names, inputfile)

    if plot_pvalues:
        compare_marginals(data_real, data_netw, t_real, t_netw, names, inputfile)

    if plot_umap:
        # Remove stimulus
        data_real = np.delete(data_real, 1, axis=0)
        data_netw = np.delete(data_netw, 1, axis=0)
        # Remove Sparc gene (index = 36)
        if p == "tests/Semrau/":
            data_real = np.delete(data_real, 36, axis=0)
            data_netw = np.delete(data_netw, 36, axis=0)
        plot_data_umap(data_real, data_netw, t_real, t_netw, inputfile)



if __name__ == "__main__":
   main(sys.argv[1:])
