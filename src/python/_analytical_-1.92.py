# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/_output/"
index = 2


# a = 1.05


target = LOCATION + "SingleAve_-1.913000_0.040000_0.140000.txt"
N = -1.913
S_min = 0.04
S_max = 0.14
n = 181
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)
# a, e = NStoAE(data[:, 0], data[:, 1])


m = data.shape[0] / n
print n, m
x = np.linspace(-180, 180, n)
x = np.linspace(-180, 180, n) / 180.0 * np.pi
y = np.linspace(S_min, S_max, m)
# e = NStoAE(N, y)[1]
X, Y = np.meshgrid(x, y)
Z = data
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)
Zmax = np.max(Z)
Zmin = np.min(Z)
level = np.linspace(Zmin, Zmin + 0.00125, 6)
level = np.insert(level, 1, Zmin + 0.000005)
# level = np.hstack((level, Zmin + 0.00047))

level2 = np.linspace(Zmin, Zmin + 0.008, 50)
# level2 = np.hstack((level2, Zmin + 0.0037))

# fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_ylim(0.3, 0.5)
# fig2, ax2 = plt.subplots()
ax.contour(X, Y, Z, level, colors='black', linewidths=1, alpha=0.25)
cmap = plt.get_cmap('magma_r')
new_cmap = truncate_colormap(cmap, 0, 0.8)
ax.contourf(X, Y, Z, level2, cmap=new_cmap)
ax.grid(color='black', linestyle='dotted', alpha=0.2, linewidth=1)
ax.grid(True)
# ax.yaxis.set_ticks(np.arange(0.05, 0.25, 0.05))
# ax.xaxis.grid(True)

# ax.contour(X, Y, Z, level2, colors='k')

# cbar = plt.colorbar(CS)x
# fig.savefig('N=-2.02.png', dpi=500, transparent=True)

plt.show()
