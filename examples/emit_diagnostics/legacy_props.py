#!/usr/bin/env python

import sys, os
import numpy as np
import h5py

def main():
    pfile = sys.argv[1]
    print('reading file ', pfile)

    h5 = h5py.File(pfile, 'r')
    print('keys: ', h5.keys())
    mass = h5.get('mass')[()]
    print('mass: ', mass)
    pz = h5.get('pz')[()]
    print('pz: ', pz)
    betagamma = pz/mass
    gamma = np.sqrt(betagamma**1 + 1)
    beta = betagamma/gamma
    print('beta: ', beta)
    print('gamma: ', gamma)

if __name__ == "__main__":
    main()
