# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes


def readData(data):
    return [data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]]


def calI(H, a, e):
    return np.arccos(H / np.sqrt(a * (1 - e**2))) * 180 / np.pi

# LOCATION = "assets/_output/"
# target = LOCATION + "print.txt"
# data = np.loadtxt(target)

# [t, x, y, z, dist] = readData(data)
# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# ax1.plot(x, y)
# ax2.plot(t, dist)

# plt.show()

print calI(-0.95, 1.00, 0.24)
