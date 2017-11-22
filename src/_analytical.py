# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/_output/"
index = 2

N = -1.98


target = LOCATION + "DoubleAve_-0.870000.txt"
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)
# x = np.arange(-180, 181, 1) / 180.0 * np.pi
n = 361
m = data.shape[0] / n
x = np.linspace(0, 360, n)
y = np.linspace(0, 0.489, m)
X, Y = np.meshgrid(x, y)
# Z = data[:, 0] - 1 * data[:, 1]
Z = data
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)

level = np.linspace(np.min(Z), np.min(Z) + 0.00013, 15)
level2 = [0.00093911]
print level
print level2
# level = np.linspace(0, 0.01, 40)

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
fig, ax = plt.subplots()

ax.contourf(X, Y, Z, level, cmap="viridis_r")
CS = ax.contour(X, Y, Z, level2, cmap="viridis_r")
# ax.plot(x, Z[489, :])
# ax.set_xlim(-180, 180)
# ax.set_ylim(0, 0.2)
# cbar = plt.colorbar(CS)x

plt.show()
