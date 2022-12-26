import sys; sys.path += ['../']
import numpy as np
from cardamom import NetworkModel
import getopt

verb = 1

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:", ["ifile="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            inputfile = arg

    p = '{}/'.format(inputfile)  # Name of the file where are the data

    print('infer {}'.format(inputfile))
    fname = p+'Data/panel_real.txt'.format(inputfile)
    data = np.loadtxt(fname, dtype=int, delimiter='\t')[1:, 1:].T
    data[:, 0] = np.loadtxt(fname, dtype=int, delimiter='\t')[0, 1:]
    G = np.size(data, 1)
    model = NetworkModel(G-1)
    model.fit(data, verb=verb)

    np.save(p+'cardamom/basal', model.basal)
    np.save(p + 'cardamom/inter', model.inter)
    np.save(p+'cardamom/basal_t', model.basal_t)
    np.save(p+'cardamom/inter_t', model.inter_t)
    np.save(p+'cardamom/kmin', np.min(model.data_bool, 0))
    np.save(p+'cardamom/kmax', np.max(model.data_bool, 0))
    np.save(p+'cardamom/data_bool', model.data_bool)
    np.save(p+'cardamom/bet', model.a[-1, :])


if __name__ == "__main__":
   main(sys.argv[1:])

