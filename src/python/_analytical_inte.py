# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy.integrate import odeint, ode


LOCATION = "assets/_output/"
index = 2


target = LOCATION + "SingleAve_-1.913000_0.040000_0.140000.txt"
N = -1.913
S_min = 0.04
S_max = 0.14
n = 721
data = np.loadtxt(target)
num = data.shape[0]
m = num / n
print n, m
xx = np.linspace(-180, 180, n)
x_interv = xx[1] - xx[0]
yy = np.linspace(S_min, S_max, m)
y_interv = yy[1] - yy[0]
X, Y = np.meshgrid(xx, yy)
Z = np.nan_to_num(data.reshape(m, n))
delta = 1e-4


def xy2xnyn(x, y):
    if x < -180 or x > 180:
        return (0, 0)
    if y < S_min or y > S_max:
        return (0, 0)
    xn = int(round(x / x_interv)) + n / 2
    yn = int(round((y - S_min) / y_interv))
    return xn, yn


def getH(x, y):
    xn, yn = xy2xnyn(x, y)

    # layer = 1
    # x_point = []
    # y_point = []
    # z_point = []

    # for i in np.arange(-layer, layer + 1):
    #     x_point.append(xx[xn + i])
    #     y_point.append(yy[yn + i])
    #     zz_point = []
    #     for j in np.arange(-layer, layer + 1):
    #         zz_point.append(Z[xn + i, yn + j])
    #     z_point.append(zz_point)
    # f = interpolate.interp2d(x_point, y_point, z_point, kind='linear')
    return Z[xn, yn]


def H_diff(x, y):
    xn, yn = xy2xnyn(x, y)
    Hx = (Z[yn, xn + 1] - Z[yn, xn - 1]) / (2 * x_interv)
    Hy = (Z[yn + 1, xn] - Z[yn - 1, xn]) / (2 * y_interv)
    return (-Hx, -Hy)


def hamiltonian(y, t):

    x1, x2 = y[0], y[1]
    Hx, Hy = H_diff(x1, x2)
    dx1 = Hy
    dx2 = -Hx
    return [dx1, dx2]

# x = np.linspace(-150, 150, 100)
# y = []
# for i in x:
#     y.append(H_diff(i, 0.11))
t = np.linspace(0, 200000, 10000)
y0 = np.array([-111.4048, 0.07480166])
y = odeint(hamiltonian, y0, t)

print "DONE!"
co = 52.2333333333
fig, ax = plt.subplots()
ax.plot(t, y[:, 0], 'r')
ax.set_ylim(-180, 180)
ax.set_xlim(0, co * 3000)
# plt.plot(x, y)
# print Z.shape


# print H_diff(178, 0.11)
# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# ax.set_ylim(0.3, 0.5)
# fig, ax = plt.subplots()
# Zmin = np.min(Z)
# level = np.linspace(Zmin, Zmin + 0.00125, 6)
# ax.contourf(X, Y, Z, level, colors='black', linewidths=1, alpha=0.25)
fig.savefig("BZ509_ana.pdf", dpi=300)
plt.show()
