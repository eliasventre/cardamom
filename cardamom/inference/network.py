"""
Core functions for network inference.

Considered sizes
----------------
    C : number of cells
    G : number of genes including stimulus

Parameters
----------
    x : array (C, G)
        Observed mRNA levels (column 0 = time points).
    inter : dictionnary of arrays (G, G)
        Given key t, inter[t][i,j] denotes interaction i -> j at time t.
    basal : array (G,)
        Basal activity for each gene.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from numba import njit

sl = 5e-3 # pseudo l1 coefficient
p = .4

@njit
def penalization_l1(x, s):
    return (x-s/2)*(x>s) - (x+s/2)*(-x>s) + ((x**2)/(2*s))*(x<=s and -x<=s)

@njit
def grad_penalization_l1(x, s):
    return 1*(x>s) - 1*(-x>s) + (x/s)*(x<=s and -x<=s)

@njit
def penalization(Q, X, X_init, time_init, l, sc, G, cnt_move, j):

    l_inter, l_diag = l[0], l[1]
    if not time_init: Q += penalization_l1(X[-1] - X_init[j, -1], sc) * l_inter * p
    if time_init == 1: Q += penalization_l1(X[0] - X_init[j, 0], sc) * l_inter * p
    tmp_list = list(range(1, j)) + list(range(j + 1, G))
    for i in tmp_list:
        tmp = cnt_move[i] * (1 + abs(X_init[i, j]) / (1 + abs(X_init[j, i])))
        Q += penalization_l1(X[i] - X_init[j, i], sc) * l_inter * tmp
    tmp_diag = cnt_move[j] / (1 + np.sum(np.abs(X_init[j, :-1])) - abs(X_init[j, j]))
    Q += l_diag * (penalization_l1(X[j] - X_init[j, j], sc) + (X[j] - X_init[j, j]) ** 2) * tmp_diag
    return Q

@njit
def grad_penalization(dq, X, X_init, time_init, l, sc, G, cnt_move, j):

    l_inter, l_diag = l[0], l[1]
    if not time_init: dq[-1] += grad_penalization_l1(X[-1] - X_init[j, -1], sc) * l_inter * p
    else: dq[-1] = 0
    if time_init == 1: dq[0] += grad_penalization_l1(X[0] - X_init[j, 0], sc) * l_inter * p
    else:  dq[0] = 0
    tmp_list = list(range(1, j)) + list(range(j + 1, G))
    for i in tmp_list:
        tmp = cnt_move[i] * (1 + abs(X_init[i, j]) / (1 + abs(X_init[j, i])))
        dq[i] += grad_penalization_l1(X[i] - X_init[j, i], sc) * l_inter * tmp
    tmp_diag = cnt_move[j] / (1 + np.sum(np.abs(X_init[j, :-1])) - abs(X_init[j, j]))
    dq[j] += l_diag * (grad_penalization_l1(X[j] - X_init[j, j], sc) + 2 * (X[j] - X_init[j, j])) * tmp_diag
    return dq

def objective(X, y, vect_kon, ko, X_init, time_init, l, sc, G, cnt_move, j):
    """
    Objective function to be maximized (all cells, one gene, all timepoints).
    """
    sigma = expit(X[-1] + y @ X[:-1])
    Q = np.sum((ko + (1 - ko) * sigma - vect_kon) ** 2)
    return penalization(Q, X, X_init, time_init, l, sc, G, cnt_move, j)


def grad_theta(X, y, vect_kon, ko, X_init, time_init, l, sc, G, cnt_move, j):
    """
    Objective gradient for gene i for all cells.
    """

    dq = np.zeros(G + 1)
    sigma = expit(X[-1] + y @ X[:-1])
    tmp = 2 * (ko + (1 - ko) * sigma - vect_kon) * (1 - ko) * sigma * (1 - sigma)
    for i in range(0, G):
        dq[i] += np.sum(y[:, i] * tmp)
    dq[-1] += np.sum(tmp)

    return grad_penalization(dq, X, X_init, time_init, l, sc, G, cnt_move, j)


def build_cnt(cnt, cnt_move, vect_kon, vect_t, times, G):

    for j in range(1, G):
        for i in range(1, G):
            tmp_cnt = 0
            if cnt > 2: # -2 for letting the first timepoints after the stimulus unchanged
                for t in range(0, cnt - 2):
                    tmp = np.sum(vect_t == times[t])
                    tmp1 = np.sum(vect_kon[vect_t == times[t], i] == 1) / tmp
                    tmp2 = np.sum(vect_kon[vect_t == times[t], j] == 1) / tmp
                    if tmp1 > p:
                        tmp_cnt += tmp1 * tmp2
                tmp_cnt /= cnt-1
            cnt_move[j, i] = (1 + tmp_cnt) * max(p, min(1 +
                np.sum(vect_kon[vect_t == times[cnt - 2], i] == 1) / np.sum(vect_t == times[cnt - 2])
              - np.sum(vect_kon[vect_t == times[cnt - 1], i] == 1) / np.sum(vect_t == times[cnt - 1]), 1))
    return cnt_move


def core_inference(variations, times, vect_t, y, vect_kon, ko, Xo, l, sl, G, cnt_init, cnt_end):

    cnt_move = np.ones((G, G))
    theta_t = np.zeros((cnt_end - cnt_init, G, G + 1))
    for cnt, time in enumerate(times[cnt_init:cnt_end]):
        p = l * np.sum(vect_t == time) * 2.5 / 100
        if cnt + cnt_init > 1:
            cnt_move = build_cnt(cnt + cnt_init, cnt_move, vect_kon, vect_t, times, G)

        X_init = Xo.copy()
        for j in range(1, G):
            res = minimize(objective, X_init[j, :],
                               args=(y[vect_t == time, :],
                                     vect_kon[vect_t == time, j],
                                     ko[j], X_init, cnt + cnt_init, p, sl, G, cnt_move[j], j),
                               jac=grad_theta,
                               method='L-BFGS-B',
                               tol=1e-5)

            if not res.success: print('Warning, minimization time {} failed'.format(time))
            variations[cnt][j, :] += np.abs(res.x[:-1] - Xo[j, :-1])
            Xo[j, :] = res.x[:]
            theta_t[cnt, j, :] = Xo[j, :]

    return Xo, variations, theta_t


def inference_optim(vect_t, times, y, vect_kon, ko):
    """
    Network inference procedure.
    Return the inferred network (basal + network) and the time at which each edge has been detected with strongest intensity.
    """
    G = np.size(y, 1)
    inter = np.zeros((G, G))
    basal = np.zeros(G)
    inter_t = np.zeros((len(times), G, G))
    basal_t = np.zeros((len(times), G))
    variations = np.zeros((len(times), G, G))
    time_variations = np.zeros((G, G))

    # penalites
    l = np.array([1, np.sqrt(G - 1)])  # l_inter, l_basal, l_diag

    Xo = np.zeros((G, G + 1))
    Xo, variations, theta_t = core_inference(variations, times, vect_t, y, vect_kon, ko, Xo, l, sl, G, 0, len(times))

    for i in range(0, G):
        for j in range(0, G):
            time_variations[i, j] = np.argmax(variations[1:, i, j])

    inter[:, :] = Xo[:, :-1]
    basal[:] = Xo[:, -1]
    inter_t[:, :, :] = theta_t[:, :, :-1]
    basal_t[:, :] = theta_t[:, :, -1]

    return basal, inter, time_variations, basal_t, inter_t