# -*- coding: UTF-8 -*-
import numpy as np


def ifunc(C, e):
    result = (C - 1 / a) / (2 * np.sqrt(a * (1 - e**2)))
    # if (result > 1.0 or result < -1.0):
    #     return np.nan
    return np.arccos(result) / np.pi * 180


a = 1
nx, ny = (101, 97)
C = np.linspace(-1.0, 1.0, nx)
E = np.linspace(0.02, 0.98, ny)

egrid, cgrid = np.meshgrid(E, C)
igrid = ifunc(cgrid, egrid)
print cgrid, egrid
cgrid = cgrid.flatten()
egrid = egrid.flatten()
igrid = igrid.flatten()


data = np.empty((0, 7))
for (c, e, i) in zip(cgrid, egrid, igrid):
    if (not np.isnan(i) and i > 90):
        data = np.vstack((data, np.array([[a, e, i, 0, 0, 0, c]])))

np.savetxt("assets/input.in", data, fmt='%.2f %.4f %10.6f %d %d %d %5.2f')
