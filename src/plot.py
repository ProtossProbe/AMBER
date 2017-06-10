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


target = ["../assets/rot.txt", "../assets/init.txt",
          "../assets/elements.txt", "../assets/key.txt"]

# data2 = np.loadtxt(target[2])
# plt.plot(data2[:, 0], data2[:, 1])
# plt.show()

data3 = np.loadtxt(target[3])
data3[:, 3] = np.unwrap(data3[:, 3])
ax = plt.subplot(111)

ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)

ax.scatter(15.874 * data3[:, 2] * np.cos(data3[:, 3]),
           15.874 * data3[:, 2] * np.sin(data3[:, 3]), s=0.5)
# plt.plot(data3[:, 0], data3[:, 1])
ax.grid(linestyle='dashed')
plt.show()
