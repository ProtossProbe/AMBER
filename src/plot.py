# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mu = 0.001


def plot_primaries(ax):
    ax.plot([1 - mu], [0], 'oy')
    ax.plot([-mu], [0], 'ob')


def orbit_plot(y, range=0.5):
    """
    给定一系列的轨道状态，画出三维空间下的轨道图
    """
    # 初始化2D绘图模块
    fig, ax = plt.subplots()

    # 指定绘图范围
    ax.set_xlim(1 - range, 1 + range)
    ax.set_ylim(-range, range)

    # 画出轨道
    if type(y) is np.ndarray:
        ax.plot(y[:, 1], y[:, 2], alpha=0.7, linewidth=1.5)
    elif type(y) is list:
        for i in y:
            ax.plot(i[:, 1], i[:, 2], alpha=0.7, linewidth=1.5)

    plot_primaries(ax)
    plt.show()


def contour_plot(target):
    fig, ax = plt.subplots()
    # 指定绘图范围
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    for file in target:
        data = np.loadtxt(file)
        ax.scatter(data[:, 2] * np.cos(data[:, 3] / 2),
                   data[:, 2] * np.sin(data[:, 3] / 2), s=3)
    ax.grid(linestyle='dashed')
    plt.show()


# target = ["../assets/rot.txt", "../assets/init.txt",
#           "../assets/elements.txt", "../assets/key.txt", "../assets/key3.txt"]

target = ["assets/out.txt"]

# target2 = ["../assets/170629/-1.03 | 0.06.txt",
#            "../assets/170629/-1.03 | 0.06-.txt",
#            "../assets/170629/-1.03 | 0.15.txt",
#            "../assets/170629/-1.03 | 0.15-.txt",
#            "../assets/170629/-1.03 | 0.2.txt",
#            "../assets/170629/-1.03 | 0.03-.txt",
#            "../assets/170629/-1.03 | 0.1-.txt",
#            "../assets/170629/-1.03 | 0.08.txt",
#            "../assets/170629/-1.03 | 0.055.txt"]


data1 = np.loadtxt(target[0])
# plt.plot(data2[:, 0], data2[:, 1])
# plt.show()
# data1 = np.loadtxt(target[2])
# plt.plot(data1[:1000, 1], data1[:1000, 2])
fig1, ax1 = plt.subplots()
fig2, [ax3, ax4] = plt.subplots(2, sharex=True)

ax1.set_xlim(-1.5, 1.5)
ax1.set_ylim(-1.5, 1.5)

# for file in target:
ax1.plot(data1[:, 1], data1[:, 2], linewidth=0.1)

ax3.plot(data1[:, 0], data1[:, 5])
ax4.plot(data1[:, 0], data1[:, 6])
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
ax1.grid(linestyle='dashed')
ax3.grid(linestyle='dashed')
ax4.grid(linestyle='dashed')
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
