"""Main class for network inference and simulation."""
import numpy as np
from cardamom.inference import (inference_optim, infer_kinetics,
    log_gamma_poisson_pdf)

np.set_printoptions(formatter={'float': '{0:.3f}'.format})


class NetworkModel:
    """Handle networks within the package."""

    def __init__(self, n_genes=None):
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
            n = n_genes + 1  # Genes plus stimulus
            # Default degradation rates
            self.d = np.zeros((2, n))
            self.d[0] = np.log(2)/9  # mRNA degradation rates
            self.d[1] = np.log(2)/46  # protein degradation rates
            # Default network parameters
            self.basal = np.zeros(n)
            self.inter = np.zeros((n, n))
            self.variations = np.zeros((n, n))

    def core_basins_binary(self, data, alph, b, g):
        """Compute the basal parameters of filtered genes."""
        for n, c in enumerate(data):
            for z, k in enumerate(alph):
                self.weight[n, g, z] = log_gamma_poisson_pdf(c, k, b)
            self.data_bool[n, g] = alph[np.argmax(self.weight[n, g, :])]

    def core_optim(self, x, vect_t, times, n_tot):
        """
        Fit the network model to the data.

        Return the list of successive objective function values.
        """
        # Initialization
        k_max = np.max(x, 0)
        ko = np.min(x, 0) / k_max
        y, vect_kon = x / k_max, x / k_max

        variations_tot = np.zeros((n_tot, n_tot))
        basal_tot = np.zeros(n_tot)
        inter_tot = np.zeros((n_tot, n_tot))
        basal_tot_t = np.zeros((len(times), n_tot))
        inter_tot_t = np.zeros((len(times), n_tot, n_tot))

        # Inference procedure
        basal, inter, variations_time, basal_t, inter_t = inference_optim(
            vect_t, times, y, vect_kon, ko)

        # Build the results
        cnt_i = 0
        for i in range(0, n_tot):
            cnt_j = 0
            if self.mask[i]:
                for j in range(0, n_tot):
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
        n_cells, g_tot = data.shape
        vect_t = data[:, 0]
        times = list(set(vect_t))
        times.sort()

        # Get kinetic parameters
        seuil = 1e-3
        self.data_bool = np.ones_like(data, dtype='float')
        self.data_bool[vect_t == 0, 0] = 0
        self.weight = np.zeros((n_cells, g_tot, 2), dtype='float')
        a = np.ones((3, g_tot))
        for g in range(1, g_tot):
            x = data[:, g]
            at, a[-1, g] = infer_kinetics(x, vect_t, verb=verb)
            a[0, g] = max(np.min(at), seuil)
            a[1, g] = max(np.max(at), seuil)
            if verb:
                print(f'Gene {g} calibrated...', a[:, g])
            self.core_basins_binary(x, a[:-1, g], a[-1, g], g)
        self.a = a

        # Remove genes with too small variations
        threshold = 0.1
        self.mask = np.ones(g_tot, dtype='bool')
        for g in range(1, g_tot):
            meang = [np.mean(data[vect_t == time, g]) for time in times]
            if np.max(meang) - np.min(meang) < threshold:
                self.mask[g] = 0

        if verb:
            print('number genes of interest', np.sum(self.mask))
        self.d = self.d[:, self.mask]
        data_bool = self.data_bool[:, self.mask]

        self.basal, self.inter, self.variations, self.basal_t, self.inter_t = \
            self.core_optim(data_bool, data[:, 0], times, g_tot)

        if verb:
            print('TOT', self.inter.T, self.basal)
