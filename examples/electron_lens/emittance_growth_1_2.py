#!/usr/bin/env python
import sys,os
import numpy as np
import tables
import matplotlib.pyplot as plt

if __name__ == "__main__":
    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    file1 = "/data/egstern/elens/elens_1.0e11.00/"
    file2 = "/data/egstern/elens/elens_2.0e11.00/"

    h5a = tables.openFile(file1 + "cell_full2_0.h5")
    h5b = tables.openFile(file2 + "cell_full2_0.h5")

    emitxa = h5a.root.emitx.read()
    emitya = h5a.root.emity.read()
    emitxb = h5b.root.emitx.read()
    emityb = h5b.root.emity.read()

    nsteps = int(emitxa.shape[0]/12.)
    xsteps = np.arange(emitxa.shape[0],dtype='d')/12.0

    plt.figure()
    plt.subplot(211)
    plt.title("x emittance")
    plt.plot(xsteps, emitxa, label='xi ~ 0.5', lw=2)
    plt.plot(xsteps, emitxb, label='xi ~ 0.9', lw=2)
    plt.xlabel("turn")
    plt.ylabel('x emittance')
    plt.legend(loc='best')

    plt.subplot(212)
    plt.title("y emittance")
    plt.plot(xsteps, emitya, label='xi ~ 0.5', lw=2)
    plt.plot(xsteps, emityb, label='xi ~ 0.9', lw=2)
    plt.xlabel("turn")
    plt.ylabel("y emittance")
    plt.legend(loc='best')

    plt.show()
