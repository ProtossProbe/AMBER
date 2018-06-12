# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

LOCATION = "assets/mercury/"


def wrapTo180(phi):
    phi = phi % 360
    for i in np.arange(phi.size):
        if phi[i] > 180:
            phi[i] -= 360
    return phi


def wrapTo360(phi):
    return phi % 360


def readElements(data, left, right):
    return [data[left:right, 0], data[left:right, 4], data[left:right, 5], data[left:right, 6], data[left:right, 7],
            data[left:right, 8], data[left:right, 9]]


def calN(a, e, I):
    I = I / 180.0 * np.pi
    return np.sqrt(a) * (np.sqrt(1 - e**2) * np.cos(I) - 1)


def calS(a, e):
    return np.sqrt(a) * (1 - np.sqrt(1 - e * e))


def poinSect_M(data_a, data_p):
    row, col = data_a.shape
    result_a = np.empty(shape=[0, col])
    result_p = result_a
    focus = data_a[:, 9]
    for i in np.arange(row - 1):
        if focus[i + 1] - focus[i] < 0:
            if focus[i + 1] < np.abs(360 - focus[i]):
                result_a = np.vstack((result_a, data_a[i + 1]))
                result_p = np.vstack((result_p, data_p[i + 1]))
            else:
                result_a = np.vstack((result_a, data_a[i]))
                result_p = np.vstack((result_p, data_p[i]))
    return result_a, result_p


mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 0.5

# mu = 0.001
target1 = ["SV131.txt"]
target = LOCATION + target1[0]
eMax = 0.819152
data_ast = np.loadtxt(target, skiprows=4)
num = data_ast.shape[0]
print num
[t1, a1, e1, I1, ome1, Ome1, M1] = readElements(
    data_ast, 0, num)
fig1, ax1 = plt.subplots()
ome1 = wrapTo180(ome1)
H = np.sqrt(a1*(1-e1**2)) * np.cos(I1/180*np.pi)
ax1.scatter(ome1, e1, s=1)
# ax1.scatter(ome1, I1)
ax1.set_xlim(-180, 180)
ax1.set_ylim(0, eMax)
# ax1.scatter(t1, e1)
# ax1.scatter(t1, I1)
# target2 = "JUPITER.txt"
# target3 = "SATURN.txt"
# data_pla = np.loadtxt(target2, skiprows=4)
# end_t = data_pla[-1, 0]
# for i in np.arange(0, 1):
#     print i
#     target = target1[i]
#     data_ast = np.loadtxt(target, skiprows=4)
#     num = data_ast.shape[0]
#     print num
#     [t1, a1, e1, I1, ome1, Ome1, M1] = readElements(
#         data_ast, 0, num)

#     [t2, a2, e2, I2, ome2, Ome2, M2] = readElements(
#         data_pla, 0, num)
#     N = calN(a1 / a2, e1, 180)
#     S = calS(a1 / a2, e1)
#     num = N.shape[0]
#     print N[num / 2]
#     t1 = t1
#     lambda_ast = ome1 + Ome1 + M1
#     lambda_pla = ome2 + Ome2 + M2

#     phi1 = wrapTo180(lambda_ast - lambda_pla - 2 * ome1)
#     phi2 = wrapTo180(lambda_ast - lambda_pla)
#     omega = wrapTo180(ome1)

#     # fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'))
#     fig1, ax1 = plt.subplots()
#     fig2, (ax2, ax3, ax4, ax5, ax6) = plt.subplots(5, sharex=True)

#     # ax1.scatter(phi1 / 180 * np.pi, e1, alpha=0.5)
#     # ax1.set_ylim(0.3, 0.5)\
#     ax1.scatter(t1, phi1, s=5, alpha=0.03)
#     ax1.set_ylim(-180, 180)
#     ax1.set_xlim(0, 3000)
#     # ax2.set_ylim(-180, 180)
#     # ax2.scatter(data_ast[:, 0], wrapToPi(phi1), alpha=0.5)
#     # ax2.scatter(data_ast[:, 0], wrapToPi(phi2), alpha=0.5)
#     # ax2.scatter(data_ast[:, 0], wrapToPi(phi2), alpha=0.5)
#     ax2.scatter(t1, N, alpha=0.5)
#     ax3.scatter(t1, a1, alpha=0.5)
#     ax4.scatter(t1, S, alpha=0.5)
#     ax5.scatter(t1, I1, alpha=0.5)
#     ax6.scatter(t1, omega, alpha=0.5)

#     for axx in [ax2, ax3, ax4, ax5, ax6]:
#         axx.grid(color='grey', linestyle='dashed')
#     ax1.set_axis_off()
#     fig1.savefig("BZ509.png", dpi=300, transparent=True)
#     # fig2.savefig("_A_" + str(i + 1) + "_2.png", dpi=150)
#     print a1[0] / a2[0], phi1[0], S[0], N[0]
plt.show()
