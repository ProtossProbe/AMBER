# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import _utility as ut

mu = 0.001
mu_s = 1 - mu
rH = (mu / (3 * mu_s)) ** (1. / 3.)

Jupiter_a = 5.20
Saturn_a = 9.58
Uranus_a = 19.20
Neptune_a = 30.05

jpl_data = np.genfromtxt('results_q>2.csv', delimiter=',', skip_header=1)

# [a, e, I, Ome, ome, q, Q] = jpl_data[:, 1:7]
a = jpl_data[:, 1]
e = jpl_data[:, 2]
I = jpl_data[:, 3]
Ome = jpl_data[:, 4]
ome = jpl_data[:, 5]
q = jpl_data[:, 6]
Q = jpl_data[:, 7]

fig, ax = plt.subplots()
ax.scatter(a, I, s=10)
# ax.hist(a, bins=72, range=[5, 50])
ax.set_xscale('log')
ax.plot([Jupiter_a, Jupiter_a], [0, 180], color="red")
ax.plot([Saturn_a, Saturn_a], [0, 180], color="red")
ax.plot([Uranus_a, Uranus_a], [0, 180], color="red")
ax.plot([Neptune_a, Neptune_a], [0, 180], color="red")
ax.set_ylim(90, 180)
# ax.set_ylim(0, 1)
# ax.set_xlim(5, 250)
ax.grid(linestyle="dashed")

fig, ax2 = plt.subplots()
ax2.hist(I, bins=18, range=[90, 180])
# ax2.hist(ut.wrapTo360(ome-Ome), bins=36, range=[0, 360])
ax2.grid(linestyle="dashed")


plt.show()
