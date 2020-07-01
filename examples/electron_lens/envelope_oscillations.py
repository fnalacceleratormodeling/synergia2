#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt

import tables

os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

file_nosc = "/data/egstern/elens/bunch_rotate_nosc.00/cell_full2_0.h5"
file_2e11 = "/data/egstern/elens/bunch_rotate_2.0e11.00/cell_full2_0.h5"

h5_nosc = tables.openFile(file_nosc)
h5_2e11 = tables.openFile(file_2e11)

plt.subplot(221)
plt.title("<x,x>")
plt.plot(h5_nosc.root.std[0,:]**2, 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.std[0,:]**2, 'o-', label="xi ~0.9", lw=2)
plt.xlabel("cell")
plt.legend(loc='best')

plt.subplot(222)
plt.title("<x, xp>")
plt.plot(h5_nosc.root.mom2[0,1,:], 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.mom2[0,1,:], 'o-', label="xi ~0.9", lw=2)
plt.xlabel("cell", fontsize='large')
plt.legend(loc='best')

plt.subplot(223)
plt.title("<x, y>")
plt.plot(h5_nosc.root.mom2[0,2,:], 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.mom2[0,2,:], 'o-', label="xi ~0.9", lw=2)
plt.xlabel("cell", fontsize='large')
plt.legend(loc='best')

plt.subplot(224)
plt.title("<y, yp>")
plt.plot(h5_nosc.root.mom2[2,3,:], 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.mom2[2,3,:], 'o-', label="xi ~0.9", lw=2)
plt.xlabel("cell", fontsize='large')
plt.legend(loc='best')

plt.figure()

plt.subplot(211)
plt.title("emit x")
plt.plot(h5_nosc.root.emitx, 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.emitx, 'o-', label="xi ~0.9", lw=2)
ax = plt.gca()
ax.set_ylim([0.980e-6, 1.04e-6])
plt.xlabel("cell")
plt.ylabel("x emittance")
plt.legend(loc='best')

plt.subplot(212)
plt.title("emit y")
plt.plot(h5_nosc.root.emity, 'o-', label="no SC", lw=2)
plt.plot(h5_2e11.root.emity, 'o-', label="xi ~0.9", lw=2)
ax = plt.gca()
ax.set_ylim([0.980e-6, 1.04e-6])
plt.xlabel("cell")
plt.ylabel("y emittance")
plt.legend(loc='best')

plt.show()
