# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/_output/"
index = 2

N = -1.98


target = LOCATION + "DoubleAve_-0.950000.txt"
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)

n = 361
m = data.shape[0] / n
print n, m
x = np.linspace(-180, 180, n)
x = np.linspace(-180, 180, n) / 180.0 * np.pi
y = np.linspace(0, 0.06, m)
X, Y = np.meshgrid(x, y)
Z = data
# Z = data[:, 0] - 1 * data[:, 1]
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)

level = np.linspace(np.min(Z), np.min(Z) + 0.00015, 20)
# level2 = [0.00093911]
# level = np.linspace(0, 0.01, 40)


# fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# fig2, ax2 = plt.subplots()
ax.contourf(X, Y, Z, level, cmap="viridis_r")
ax.grid(False)
# CS = ax.contour(X, Y, Z, level2, cmap="viridis_r")
# ax.plot(x, Z[489, :])
# ax.set_xlim(-180, 180)
# ax.set_ylim(0, 0.02)
# cbar = plt.colorbar(CS)x

plt.show()
