# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import _utility as ut


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/old/1206/planar/"
index = 2


# a = 1.05


NAME = "SingleAve_-1.96_0.08"
target = LOCATION + NAME + ".txt"
words = NAME.split("_")
N = float(words[1])
S_min = 0
S_max = float(words[2])
print N, S_min, S_max
n = 721
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)
# a, e = NStoAE(data[:, 0], data[:, 1])


m = data.shape[0] / n
print n, m
# x = np.linspace(-180, 180, n)
x = np.linspace(0, 360, n) / 180.0 * np.pi
y = np.linspace(S_min, S_max, m)
# y = np.sqrt(2*y)
y = ut.NStoAE(N, y)[1]
X, Y = np.meshgrid(x, y)
Z = data
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)
Zmax = np.max(Z)
Zmin = np.min(Z)
# level = np.linspace(Zmin, Zmin + 0.0002, 10)
# level = np.insert(level, 1, Zmin + 0.000005)
# level = np.hstack((level, Zmin + 0.00047))

level2 = np.linspace(Zmin, Zmin + 0.002, 40)
# level2 = np.linspace(Zmax - 0.010, Zmax, 30)
# level2 = np.hstack((level2, Zmin + 0.0037))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'))
# ax2.set_ylim(0, 0.008)
# fig2, ax2 = plt.subplots()
# ax2.contour(X, Y, Z, level, colors='black', linewidths=1, alpha=0.25)
# cmap = plt.get_cmap("gist_rainbow")
ax1.contour(X, Y, Z, level2)
ax2.contour(X, Y, Z, level2)
ax2.grid(color='black', linestyle='dotted', alpha=0.2, linewidth=1)
ax2.grid(True)
# ax.yaxis.set_ticks(np.arange(0.05, 0.25, 0.05))
# ax.xaxis.grid(True)

# ax.contour(X, Y, Z, level2, colors='k')

# cbar = plt.colorbar(CS)x
# fig.savefig('N=-2.02.png', dpi=500, transparent=True)

plt.show()
