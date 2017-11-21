# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/_output/"
index = 2

N = -1.98


target = LOCATION + "disturb_-1.92.txt"
# target = LOCATION + "Ast_" + str(index) + ".txt"
# target2 = LOCATION + "A_" + str(index) + ".txt"
# data2 = np.loadtxt(target2, skiprows=4)

# levels = np.linspace(-1.6, -1.4, 20)
# print levels

data = np.loadtxt(target)
# x = np.arange(-180, 181, 1) / 180.0 * np.pi
x = np.arange(-180, 181, 1)
y = np.arange(0, 0.1, 0.0002)
n = x.shape[0]
m = y.shape[0]
X, Y = np.meshgrid(x, y)
Z = data[:, 0] - 1 * data[:, 1]
Z = Z.reshape(m, n)

level = np.linspace(np.min(Z), np.min(Z) + 0.0020, 25)
# level = np.linspace(0, 0.01, 40)

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
fig, ax = plt.subplots()

CS = ax.contourf(X, Y, Z, level, cmap='hot')
ax.grid(False)
# ax.set_xlim(-180, 180)
ax.set_ylim(0, 0.1)
# cbar = plt.colorbar(CS)x

plt.show()