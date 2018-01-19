# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes


mu = 0.001
mu_s = 1 - mu
rH = (mu / (3 * mu_s)) ** (1. / 3.)
LOCATION = "assets/_output/"


def readPos(data, num):
    return [data[:num, 0], data[:num, 1], data[:num, 2], data[:num, 3]]

name = ["Ast_1.txt", "Ast_2.txt", "Ast_3.txt", "Ast_4.txt"]
# plt.axes().set_aspect('equal', 'datalim')
for index in range(3, 4):
    target = LOCATION + name[index - 1]
    data = np.loadtxt(target)

    fig, ax = plt.subplots()

    num = data.shape[0]
    [t, x, y, z] = readPos(data, num)
    width = 0.5
    ax.plot(x, y, linewidth=width)

    ax.plot([-mu, 1 - mu], [0, 0], 'x', color='black')
    hill = plt.Circle((1 - mu, 0), rH, color='g',
                      fill=False, linewidth=1.5)
    fig.gca().add_artist(hill)
    # ax.set_ylim(0.971, 0.986)
    # plt.plot(t, lam_p)
    # fig.patch.set_visible(False)
    # ax.axis('off')
    # fig.set_size_inches(6, 6)

    # plt.plot(data1[:, 0], data1[:, 1])
    ax.set_aspect(1)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    # ax.set_xlim(0.90, 1.1)
    # ax.set_ylim(-0.1, 0.1)
    # ax.set_ylim(-0.1, 0.1)
    # fig.savefig('frame' + str(n) + '.png', dpi=100, transparent=True)

plt.show()
