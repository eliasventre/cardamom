"""
Main class for network inference and simulation
"""
import numpy as np
from ..inference import inference_optim, infer_kinetics, log_gamma_poisson_pdf

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

class NetworkModel:
    """
    Handle networks within the package.
    """
    def __init__(self, n_genes=None, times=None):
        # Kinetic parameters
        self.d = None
        # Mixture parameters
        self.data_bool = None
        # Optim parameters
        self.n_train = None
        # Network parameters
        self.mask = None
        self.a = None
        self.basal = None
        self.inter = None
        self.basal_t = None
        self.inter_t = None
        self.variations = None
        # Default behaviour
        if n_genes is not None:
            G = n_genes + 1 # Genes plus stimulus
            # Default degradation rates
            self.d = np.zeros((2,G))
            self.d[0] = np.log(2)/9 # mRNA degradation rates
            self.d[1] = np.log(2)/46 # protein degradation rates
            # Default network parameters
            self.basal = np.zeros(G)
            self.inter = np.zeros((G,G))
            self.variations = np.zeros((G, G))


    def core_basins_binary(self, data, alph, b, g):
        """
        Compute the basal parameters of filtered genes.
        """
        for n, c in enumerate(data):
            for z, k in enumerate(alph):
                self.weight[n, g, z] = log_gamma_poisson_pdf(c, k, b)
            self.data_bool[n, g] = alph[np.argmax(self.weight[n, g, :])]


    def core_optim(self, x, vect_t, times, G_tot):
        """
        Fit the network model to the data.
        Return the list of successive objective function values.
        """

        # Initialization
        k_max = np.max(x, 0)
        ko = np.min(x, 0) / k_max
        y, vect_kon = x / k_max, x / k_max

        variations_tot = np.zeros((G_tot, G_tot))
        basal_tot = np.zeros(G_tot)
        inter_tot = np.zeros((G_tot, G_tot))
        basal_tot_t = np.zeros((len(times), G_tot))
        inter_tot_t = np.zeros((len(times), G_tot, G_tot))

        # Inference procedure
        basal, inter, variations_time, basal_t, inter_t = inference_optim(vect_t, times, y, vect_kon, ko)

        # Build the results
        cnt_i, cnt_j = 0, 0
        for i in range(0, G_tot):
            cnt_j = 0
            if self.mask[i]:
                for j in range(0, G_tot):
                    if self.mask[j]:
                        inter_tot[i, j] = inter[cnt_j, cnt_i]
                        inter_tot_t[:, i, j] = inter_t[:, cnt_j, cnt_i]
                        variations_tot[i, j] = variations_time[cnt_j, cnt_i]
                        cnt_j += 1
                basal_tot[i] = basal[cnt_i]
                basal_tot_t[:, i] = basal_t[:, cnt_i]
                cnt_i += 1

        return basal_tot, inter_tot, variations_tot, basal_tot_t, inter_tot_t


    def fit(self, data, verb=False):
        """
        Fit the network model to the data.
        Return the list of successive objective function values.
        """
        C, G_tot = data.shape
        vect_t = data[:, 0]
        times = list(set(vect_t))
        times.sort()

        # Get kinetic parameters
        seuil = 1e-3
        self.data_bool = np.ones_like(data, dtype='float')
        self.data_bool[vect_t == 0, 0] = 0
        self.weight = np.zeros((C, G_tot, 2), dtype='float')
        a = np.ones((3, G_tot))
        for g in range(1, G_tot):
            x = data[:, g]
            at, a[-1, g] = infer_kinetics(x, vect_t, verb=verb)
            a[0, g] = max(np.min(at), seuil)
            a[1, g] = max(np.max(at), seuil)
            if verb: print('Gene {} calibrated...'.format(g), a[:, g])
            self.core_basins_binary(x, a[:-1, g], a[-1, g], g)
        self.a = a
        
        # Remove genes with too small variations
        self.mask = np.ones(G_tot, dtype='bool')
        for g in range(1, G_tot):
            meang = [np.mean(data[vect_t == time, g]) for time in times]
            if np.max(meang) - np.min(meang) < .1:
                self.mask[g] = 0

        if verb: print('number genes of interest', np.sum(self.mask))
        self.d = self.d[:, self.mask]
        data_bool = self.data_bool[:, self.mask]

        self.basal, self.inter, self.variations, self.basal_t, self.inter_t = \
            self.core_optim(data_bool, data[:, 0], times, G_tot)

        if verb: print('TOT', self.inter.T, self.basal)

