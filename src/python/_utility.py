import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mu = 0.001
mu_s = 1 - mu


def calAEI(N, a, Sz):
    S = 2 * np.sqrt(a) - (Sz - N)
    return getAEI(N, S, Sz)


def getAEI(N, S, Sz):
    L = (S + Sz - N) / 2
    e = (Sz - N - S) / (S + Sz - N)
    I = Sz / (L * e) - 1
    e = np.sqrt(1 - e * e)
    I = np.arccos(I) * 180 / np.pi
    a = L**2
    return a, e, 180 - I


def NStoAE(N, S):
    B = S + N
    C = S - N
    a = (C / 2)**2
    e = np.sqrt(1 - (B / C)**2)
    return a, e


def wrapTo180(phi):
    phi = phi % 360
    for i in np.arange(phi.size):
        if phi[i] > 180:
            phi[i] -= 360
    return phi


def wrapTo360(phi):
    return phi % 360


def L_action(a):
    return np.sqrt(a)


def G_action(a, e):
    return L_action(a) * np.sqrt(1 - e * e)


def H_action(a, e, i):
    return G_action(a, e) * np.cos(i / 180 * np.pi)


def S(a, e, i):
    return L_action(a) - G_action(a, e)


def N(a, e, i, kj, k):
    return -(kj + k) * L_action(a) / k + S(a, e, i) + Sz(a, e, i)


def Sz(a, e, i):
    return G_action(a, e) + H_action(a, e, i)


# def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
#     new_cmap = colors.LinearSegmentedColormap.from_list(
#         'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
#         cmap(np.linspace(minval, maxval, n)))
#     return new_cmap
