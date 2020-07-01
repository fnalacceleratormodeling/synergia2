#!/usr/bin/env python

import sys, os
import numpy as np
import tables
import matplotlib.pyplot as plt

def calculate_alpha_beta(m):
    mxx = m[0,0]
    mxxp = m[0,1]
    mxpxp= m[1,1]

    betasq = mxx**2/(mxx * mxpxp - mxxp**2)
    beta = np.sqrt(betasq)
    alpha = -beta * mxxp/mxx
    return (alpha, beta)

def plot_alpha_beta(full2_file):
    # h5full = tables.open_file(sys.argv[1])
    # h5trks = tables.open_file(sys.argv[2])
    #h5full = tables.open_file("cell_full2_0.h5")
    h5full = tables.open_file(full2_file)
    #h5trks = tables.open_file("bulk_tracks_0.h5")
    
    nturns = h5full.root.mom2.shape[2]
    print h5full.root.mom2.shape[2], " turns worth of moments"

    alphax = np.zeros(nturns, dtype='d')
    betax = np.zeros(nturns, dtype='d')
    alphay = np.zeros(nturns, dtype='d')
    betay = np.zeros(nturns, dtype='d')
    for turn in range(nturns):
        alpha, beta = calculate_alpha_beta(h5full.root.mom2[0:2, 0:2, turn])
        alphax[turn] = alpha
        betax[turn] = beta
        alpha, beta = calculate_alpha_beta(h5full.root.mom2[2:4, 2:4, turn])
        alphay[turn] = alpha
        betay[turn] = beta
    h5full.close()

    plt.plot(betax, label='beta x')
    plt.plot(betay, label='beta y')
    plt.xlabel('turn')
    plt.ylabel('beta')
    plt.legend(loc='best')

    plt.figure()
    plt.plot(alphax, label='alpha x')
    plt.plot(alphay, label='alpha y')
    plt.xlabel('turn')
    plt.ylabel('alpha')
    plt.legend(loc='best')

    plt.show()
    sys.exit(0)
    pnum = 1
    xtune1 = np.zeros(nturns-1)
    pa = np.zeros(nturns-1)
    tune = np.zeros(nturns-1)
    for turn in range(nturns-1):
        # get x tune for track 1
        beta = betax[turn]
        alpha = alphax[turn]
        coords0 = h5trks.root.track_coords[pnum, 0:4, turn]
        z0 = complex(coords0[0], alpha*coords0[0]+beta*coords0[1])
        coords1 = h5trks.root.track_coords[pnum, 0:4, turn+1]
        beta = betax[turn+1]
        alpha = alphax[turn+1]
        z1 = complex(coords1[0], alpha*coords1[0]+beta*coords1[1])
        phase = np.log(z0/z1).imag
        if phase < 0.0:
            phase += 2.0*np.pi
        pa[turn] = phase
        tune[turn] = phase/(2.0*np.pi)

    avgwin = 10
    avgtune = np.zeros(nturns-1-avgwin)
    for i in range(nturns-1-avgwin):
        avgtune[i] = tune[i:i+avgwin].mean()

    plt.figure()
    plt.plot(pa, label='phase advance')
    plt.xlabel('turn')
    plt.ylabel('phase advance')

    plt.figure()
    plt.plot(tune, label='tune')
    plt.xlabel('turn')
    plt.ylabel('tune')

    plt.figure()
    plt.plot(avgtune, label='avg tune')
    plt.xlabel('turn')
    plt.ylabel('avg tune')

    plt.show()

if __name__ == "__main__":

    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"
