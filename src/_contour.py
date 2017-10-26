# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def ifunc(C, e):
    result = (C - 1 / a) / (2 * np.sqrt(a * (1 - e**2)))
    # if (result > 1.0 or result < -1.0):
    #     return np.nan
    return np.arccos(result) / np.pi * 180


def efunc(C, i):
    result = 1 - ((C - 1 / a) / (np.cos(i / 180 * np.pi) * 2))**2 / a
    if result.any() >= 0:
        return np.sqrt(result)
    else:
        return np.nan


def cfunc(e, i):
    return 1 / a + 2 * np.sqrt(a * (1 - e**2)) * np.cos(i / 180 * np.pi)

a = 1
target = ["assets/1026/summary.out"]
data = np.loadtxt(target[0])
c = data[data[:, 1] <= 1, :]
x, y, z = c[:, 1], c[:, 2], c[:, 4]

plt.style.use("fivethirtyeight")

fig, ax = plt.subplots()
xi = np.linspace(0, 1, 2000)
yi = np.linspace(90, 180, 7)
print yi

for i in yi:
    zi = cfunc(xi, i)
    ax.plot(zi, xi, alpha=0.5)


im = ax.tripcolor(x, y, z, 15)
fig.colorbar(im)
ax.set_xlim(-1, 1)
ax.set_ylim(0, 1)
# ax.scatter(x, y, s=1, alpha=0.2)
ax.grid(linestyle='dashed')
plt.savefig("map.pdf")
plt.show()
