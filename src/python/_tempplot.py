# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
from scipy.interpolate import spline
from scipy.interpolate import CubicSpline
from scipy.interpolate import splprep, splev


def readData(data):
    return [data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]]


def calI(H, a, e):
    return np.arccos(H / np.sqrt(a * (1 - e**2))) * 180 / np.pi


def cala(N, e):
    return (N / (1 + np.sqrt(1 - e**2)))**2


def getAEI(N, S, Sz):
    L = (S + Sz - N) / 2
    e = (Sz - N - S) / (S + Sz - N)
    I = Sz / (L * e) - 1
    e = np.sqrt(1 - e * e)
    I = np.arccos(I) * 180 / np.pi
    a = L**2

    return a, e, I


def calSz(N, a):
    return 2 * np.sqrt(a) + N


def calAEI(N, a, Sz):
    S = 2 * np.sqrt(a) - (Sz - N)
    return getAEI(N, S, Sz)


def cale(N, a):
    e = -N / np.sqrt(a)
    e = np.sqrt(1 - (1 - e)**2)
    return e


def calN(a, e, I):
    I = I / 180.0 * np.pi
    return np.sqrt(a) * (np.sqrt(1 - e**2) * np.cos(I) - 1)

jupiter_a = 5.20336301
jupiter_i = 1.30530

saturn_a = 9.53707032
saturn_i = 2.48446

print calN(5.1388582 / jupiter_a, 0.3806146, 180)
print calN(9.6769180 / saturn_a, 0.7623217, 180)
# LOCATION = "assets/1130/"
# NAME = "EqPoints.txt"
# fig, ax = plt.subplots()
# data = np.loadtxt(LOCATION + NAME)
# N = data[:, 0]
# S = data[:, 1]
# cs = CubicSpline(N, S)

# # a, e, I = getAEI(N, S, 0)
# # ax.plot(a, e)
# N = np.linspace(-2.1, -1.0, 1000)
# S = cs(N)
# a, e, I = getAEI(N, S, 0)
# ax.plot(a, e, linewidth=3)
# a = np.linspace(0.9, 1.1, 1000)

# for N in [-2.06, -2.02, -2.04, -2.0, -1.95, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3]:
#     e = cale(N, a)
#     ax.plot(a, e, color='grey', linestyle="dashed", alpha=0.5)

# N_list = np.array([-2.06, -2.04, -2.03, -2.02, -2.01, -2, -1.99, -1.98, -
#                    1.95, -1.90, -1.80, -1.60, -1.30])
# e_min = np.array([0.075, 0.054, 0.042, 0.033, 0.024, 0.017, 0.025, 0.036,
#                   0.13, 0.316, 0.526, 0.763, 0.94])
# e_max = np.array([0.14, 0.179, 0.209, 0.243, 0.280, 0.313, 0.345, 0.375,
#                   0.445, 0.535, 0.666, 0.834, 0.966])
# N_list = N_list[::-1]
# e_min = e_min[::-1]
# e_max = e_max[::-1]
# # cs = CubicSpline(N_list, e_min)

# # e_min = cs(N)
# a_min = cala(N_list, e_min)
# a_max = cala(N_list, e_max)
# ax.plot(a_min, e_min, linewidth=3)
# ax.plot(a_max, e_max, linewidth=3)

# # tck, u = splprep([a_min, e_min], s=0)
# # u = np.linspace(0, 1, 500)
# # new_points = splev(u, tck)
# # ax.plot(new_points[0], new_points[1], linewidth=3)

# # tck, u = splprep([a_max, e_max], s=0)
# # u = np.linspace(0, 1, 500)
# # new_points = splev(u, tck)
# # ax.plot(new_points[0], new_points[1], linewidth=3)

# ax.set_xlim(0.92, 1.08)
# ax.set_ylim(0, 0.95)
# plt.show()

# LOCATION = "assets/_output/"
# target = LOCATION + "print.txt"
# data = np.loadtxt(target)

# [t, x, y, z, dist] = readData(data)
# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# ax1.plot(x, y)
# ax2.plot(t, dist)

# plt.show()

# print getAEI(-1.973, 0.0346710935637, 0)
# print calSz(-2.017, 1.02)
# print calAEI(-1.94, 1.004, 0.0033)
# for e in np.arange(0.17, 0.19, 0.001):

# for delta in np.arange(0, 360, 10):
#     print cala(-1.65, 0.75), 0.75, 180, 0 + delta, 0, -180 + delta
# fig.savefig('N=-2.02.pdf', dpi=500, transparent=True)
