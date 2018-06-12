# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# import _tempplot


def readElements(data, left, right):
    return [data[left:right, 0], data[left:right, 4], data[left:right, 5], data[left:right, 6], data[left:right, 7],
            data[left:right, 8], data[left:right, 9]]


def wrapTo180(phi):
    phi = phi % 360
    for i in np.arange(phi.size):
        if phi[i] > 180:
            phi[i] -= 360
    return phi


def calAEI(N, a, Sz):
    S = 2 * np.sqrt(a) - (Sz - N)
    return getAEI(N, S, Sz)


def getAEI(N, S, Sz):
    L = (S + Sz - N) / 2
    e = (Sz - N - S) / (S + Sz - N)
    I = Sz / (L * e) - 1
    e = np.sqrt(1 - e * e)
    I = np.arccos(I) * 180 / np.pi
    a = L**2

    return a, e, 180 - I


mu = 0.001
mu_s = 1 - mu


LOCATION = "assets/_output/"
NAME = "RingAveOme_70.000000_8.000000_0.939693"
target = LOCATION + NAME + ".txt"

words = NAME.split("_")
eMax = float(words[3])


data = np.loadtxt(target)
# a, e = NStoAE(data[:, 0], data[:, 1])

n = 361
m = data.shape[0] / n
print n, m
x = np.linspace(-180, 180, n)
# x = np.linspace(-180, 180, n) / 180.0 * np.pi
y = np.linspace(0, eMax, m)
# y = calAEI(N, a, y)[2]
X, Y = np.meshgrid(x, y)
Z = data
Z[Z == np.inf] = 0
Z = Z.reshape(m, n)
Z = np.nan_to_num(Z)
Zmax = np.max(Z)
Zmin = np.min(Z)

print Zmax, Zmin

# fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
# level = np.linspace(Zmin, Zmin+0.002, 30)
fig, ax = plt.subplots()
ax.contourf(X, Y, Z, 20)

LOCATION = "assets/mercury/"
target1 = ["SV131.txt"]
target = LOCATION + target1[0]
data_ast = np.loadtxt(target, skiprows=4)
num = data_ast.shape[0]
print num
[t1, a1, e1, I1, ome1, Ome1, M1] = readElements(
    data_ast, 0, num)

ome1 = wrapTo180(ome1)
H = np.sqrt(a1*(1-e1**2)) * np.cos(I1/180*np.pi)
ax.scatter(ome1, e1, s=1, color='r')
# ax1.scatter(ome1, I1)
ax.set_xlim(-180, 180)
ax.set_ylim(0, eMax)

plt.show()
