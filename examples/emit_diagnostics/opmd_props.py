#!/usr/bin/env python

import sys, os
import numpy as np

import openpmd_api as io

def main():
    pfile = sys.argv[1]
    print('reading file ', pfile)

    s = io.Series(pfile, io.Access_Type.read_only)
    print('IO Series properties: ', dir(s))
    print('IO Series attributes: ', s.attributes)

    i = s.iterations
    N = len(i)
    print('number of iterations: ', N)

    for k in range(N):
        it = i[k]
        print('--------------------------------------------------')
        print('iteration ', k)
        print('iteration ', k, ' properties: ', dir(it))
        print('iteration ', k, 'attributes: ', it.attributes)

        print()
        print('it: ', it)
        print('it.particles properties: ', dir(it.particles))
        print('it.particles.attributes: ', it.particles.attributes)
        print('it.particles: ', it.particles)
        print('it.particles.items(): ', it.particles.items())
        c = 0
        for itm in it.particles.items():
            print('it.particles.items(): ', c, ': ', itm)
            c = c + 1

        print()
        parts = it.particles["bunch_particles"]
        print('it.particles[bunch_particles] properties: ', dir(parts))
        print('it.particles[bunch_particles].attributes: ', parts.attributes)
        print('beta_ref: ', parts.get_attribute('beta_ref'))
        print('gamma_ref: ', parts.get_attribute('gamma_ref'))
        masks = it.particles["bunch_particles_masks"]
        print('it.particles[bunch_particles_masks] properties: ', dir(masks))
        print('it.particles[bunch_particles_masks].attributes: ', masks.attributes)
        print()

    print('--------------------------------------------------')


if __name__ == "__main__":
    main()
