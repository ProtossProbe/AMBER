# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# import _tempplot


def readElements(data):
    return [data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4],
            data[:, 5], data[:, 6], data[:, 8]]


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/1206/planar/"
index = 2


# a = 1.05


target = LOCATION + "SingleAve_-1.200000_0.600000_1.000000.txt"
N = -1.2
S_min = 0.6
S_max = 1.0
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)
# a, e = NStoAE(data[:, 0], data[:, 1])

n = 721
m = data.shape[0] / n
print n, m
x = np.linspace(0, 360, n)
x = x - 180
x = x / 180.0 * np.pi
y = np.linspace(S_min, S_max, m)
e = NStoAE(N, y)[1]
print np.max(e)
# y = calAEI(N, a, y)[2]
X, Y = np.meshgrid(x, e)
Z = data
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)
Zmax = np.max(Z)
Zmin = np.min(Z)
level = np.linspace(Zmin, Zmin + 0.0009, 7)
level = np.insert(level, 1, Zmin + 0.000005)
level = np.hstack((level, Zmin + 0.0012))

level2 = np.linspace(Zmin, Zmin + 0.0036, 100)
# level2 = np.hstack((level2, Zmin + 0.0037))

# fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# ax.set_xlim(0, np.pi / 2)
ax.set_ylim(0.971, 0.986)
# fig2, ax2 = plt.subplots()
ax.contour(X, Y, Z, level, colors='black', linewidths=1, alpha=0.25)
cmap = plt.get_cmap('magma_r')
new_cmap = truncate_colormap(cmap, 0, 0.8)
ax.contourf(X, Y, Z, level2, cmap=new_cmap)
# ax.grid(color='black', linestyle='dotted', alpha=0.2, linewidth=1)
ax.grid(False)
# ax.yaxis.set_ticks(np.arange(0.05, 0.25, 0.05))
# ax.xaxis.grid(True)

# ax.contour(X, Y, Z, level2, colors='k')

# cbar = plt.colorbar(CS)x
fig.savefig('N=-1.20.png', dpi=500, transparent=True)
LOCATION = "assets/_output/"
name = ["Ast_3.txt"]
for index in range(1, 2):
    target = LOCATION + name[index - 1]

    data = np.loadtxt(target)
    [t, a, e, i, ome, Ome, M, Megno] = readElements(data)
    H_val = H_action(a, e, i)
    S_val = S(a, e, i)
    N_val = N_v(a, e, i)
    Sz_val = Sz(a, e, i)
    ome_b = ome - Ome
    lam = ome_b + M
    lam_p = t * 360

    sig1 = wrapTo360((lam - lam_p - 2 * ome_b))
    sig2 = wrapTo360((lam - lam_p + 2 * Ome))
    # sig_fast = wrapToPi(lam - lam_p)
    # ax.scatter(ome / 180 * np.pi, Sz_val, s=0.5, alpha=0.5)
    # ax.scatter(t, ome, s=0.5, alpha=0.5)
    # ax.scatter(t, ome, s=0.5, alpha=0.5)
    # ax2.scatter(t, ome, s=0.5, alpha=0.5)
    # ax2.scatter(t, sig2, s=0.5, alpha=0.5)
    # ax2.scatter(t, sig2, s=0.5, alpha=0.5)
    # ax.scatter(t, ome, s=0.5, alpha=0.5)
    # ax.scatter(t, i, s=0.5, alpha=0.5)
    # ax.scatter(t, sig / 180 * np.pi, s=0.5, alpha=0.5)
    # ax.grid(linestyle='dashed')
    # ax.set_xlim(0, 360)
    # ax.set_ylim(0, 0.00001)
    # ax.scatter(np.cos(sig1 / 180.0 * np.pi) * e,
    #            np.sin(sig1 / 180.0 * np.pi) * e, s=0.5)
    # ax.scatter(sig1 / 180.0 * np.pi, e, s=0.5, alpha=0.01, color="green")

    print np.mean(a), np.mean(N_val), np.mean(Sz_val), ome[0], sig1[0], sig2[0]
plt.show()
