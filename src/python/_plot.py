# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes


mu = 0.001
mu_s = 1 - mu


def readElements(data):
    return [data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4],
            data[:, 5], data[:, 6], data[:, 8]]

LOCATION = "assets/0102/"

# fig, (ax, ax2, ax3) = plt.subplots(3, sharex=True)
fig, ax = plt.subplots()
# fig2, ax2 = plt.subplots()
# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# fig3, ax3 = plt.subplots()
# a = np.linspace(0.95, 1.05, 100)
# e = np.linspace(0, 0.2, 100)

# X, Y = np.meshgrid(a, e)
# Z = S(X, Y)

# level = np.linspace(np.min(Z), np.max(Z), 10)
# CS = ax.contourf(X, Y, Z, level)
# cbar = plt.colorbar(CS)

# print S(1, 0.01)
# LOCATION = "assets/1227/"
name = ["Ast_-1.20.txt", "Ast_-1.97.txt", "Ast_-2.00.txt"]
# plt.axes().set_aspect('equal', 'datalim')
for index in range(1, 2):
    target = LOCATION + name[index - 1]
    # target = LOCATION + "Ast_19_N=-1.96.txt"
    data = np.loadtxt(target)
    [t, a, e, i, ome, Ome, M, Megno] = readElements(data)
    H_val = H_action(a, e, i)
    N_val = N(a, e, i)
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
    # ax.scatter(t, e, s=0.5, alpha=0.5)
    width = 1.0
    ax.plot(t, a, linewidth=width)
    ax2.plot(t, e, linewidth=width)
    ax3.plot(t, sig1, linewidth=width)
    # ax4.plot(t, ome, linewidth=width)
    ax.plot((0, 150), (1, 1), linewidth=1.5, color="red", linestyle='dashed')
    ax3.plot((0, 150), (180, 180), linewidth=1.5,
             color="red", linestyle='dashed')
    # ax2.scatter(t, ome, s=0.5, alpha=0.5)
    # ax2.scatter(t, sig2, s=0.5, alpha=0.5)
    # ax2.scatter(t, Sz_val, s=0.5, alpha=0.5)
    # ax.scatter(t, ome, s=0.5, alpha=0.5)
    # ax.scatter(t, i, s=0.5, alpha=0.5)
    # ax.scatter(t, sig / 180 * np.pi, s=0.5, alpha=0.5)
    # ax.grid(linestyle='dashed')
    # ax.set_xlim(0, 360)
    # ax.set_ylim(0, 0.00001)
    # ax.scatter(np.cos(sig1 / 180.0 * np.pi) * e,
    #            np.sin(sig1 / 180.0 * np.pi) * e, s=0.5)
    # ax.scatter(sig1 / 180.0 * np.pi, e, s=1, alpha=0.05)

    print np.mean(a), np.mean(N_val), np.mean(Sz_val), ome[0], sig1[0], sig2[0]
xlim = 20
ax.set_xlim(0, xlim)
ax2.set_xlim(0, xlim)
ax3.set_xlim(0, xlim)
ax3.set_ylim(135, 225)

# ax.set_ylim(0.971, 0.986)
# plt.plot(t, lam_p)
# fig.patch.set_visible(False)
# ax.axis('off')
fig.set_size_inches(6, 6)
fig.savefig('numerical.pdf', dpi=500, transparent=True)
# plt.plot(data1[:, 0], data1[:, 1])


# fig1, [ax1, ax2, ax3] = plt.subplots(3, sharex=True)

# ax1.plot(data1[:, 0], data1[:, 1])
# ax1.plot(data2[:, 0], data2[:, 4])
# ax1.grid(linestyle='dashed')

# ax2.plot(data1[:, 0], data1[:, 2])
# ax2.plot(data2[:, 0], data2[:, 5])
# ax2.grid(linestyle='dashed')

# ax3.plot(data1[:, 0], data1[:, 3])
# ax3.plot(data2[:, 0], data2[:, 6])
# ax3.grid(linestyle='dashed')


# plt.plot(data2[:, 0], data2[:, 1])
# plt.show()
# data1 = np.loadtxt(target[2])
# plt.plot(data1[:1000, 1], data1[:1000, 2])
# fig1, [ax1, ax2] = plt.subplots(2)
# fig2, [ax3, ax4] = plt.subplots(2, sharex=True)

# for file in target:
# ax1.plot(data1[:, 0], data1[:, 1], linewidth=1)
# ax2.plot(data1[:, 0], data1[:, 2], linewidth=1)

# ax3.plot(data1[:, 0], data1[:, 5])
# ax1.plot(data1[:n, 1], data1[:n, 2], 'r', linewidth=2)
# plot_primaries(ax1)
# ax2.plot(data2[:, 1], data2[:, 2], linewidth=0.5)
# ax3.scatter(data4[:, 2] * np.cos(data4[:, 3]),
# data4[:, 2] * np.sin(data4[:, 3]), s=1)
# ax2.plot(data5[:, 0], data5[:, 3])
# ax2.plot(data4[:, 0], data4[:, 3])

# ax2.plot(data2[:, 0], np.unwrap(data2[:, 3] -
#                                 data2[:, 4] + data2[:, 0]), linewidth=1)

# ax3.scatter(data5[:, 2] * np.cos(data5[:, 3]),
#             data5[:, 2] * np.sin(data5[:, 3]), s=0.2)
# ax3.scatter(data4[:, 2] * np.cos(data4[:, 3]),
#             data4[:, 2] * np.sin(data4[:, 3]), s=0.5)

# ax1.scatter(data1[:, 4] * 180 / np.pi, data1[:, 1], s=0.5)
# ax2.scatter(data2[:, 2] * np.cos(data2[:, 3]),
#             data2[:, 2] * np.sin(data2[:, 3]), s=0.5)
# ax1.grid(linestyle='dashed')
# ax2.grid(linestyle='dashed')
# ax3.grid(linestyle='dashed')
# ax4.grid(linestyle='dashed')

# plt.show()

# data3 = np.loadtxt(target[3])
# data4 = np.loadtxt(target[4])
# # data3[:, 3] = np.unwrap(data3[:, 3])
# # data3[:, 4] = np.unwrap(data3[:, 4])
# fig1, ax = plt.subplots()
# fig2, ax2 = plt.subplots()


# ax.scatter(data3[:, 4] * 180 / np.pi, data3[:, 1], s=0.5)
# # plt.plot(data1[:, 0], data1[:, 2])


# ax2.grid(linestyle='dashed')

plt.show()
