# -*- coding: UTF-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

mu = 0.001
LOCATION = "assets/_output/"
index = 5

target = LOCATION + "Ast_" + str(index) + ".txt"
target2 = LOCATION + "A_" + str(index) + ".txt"
data2 = np.loadtxt(target2, skiprows=4)

data1 = np.loadtxt(target)
data1[:, 0] = data1[:, 0] * 1.00027
fig1, [ax1, ax2, ax3] = plt.subplots(3, sharex=True)

ax1.plot(data1[:, 0], data1[:, 1])
ax1.plot(data2[:, 0], data2[:, 4])
ax1.grid(linestyle='dashed')

ax2.plot(data1[:, 0], data1[:, 2])
ax2.plot(data2[:, 0], data2[:, 5])
ax2.grid(linestyle='dashed')

ax3.plot(data1[:, 0], data1[:, 3])
ax3.plot(data2[:, 0], data2[:, 6])
ax3.grid(linestyle='dashed')


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
