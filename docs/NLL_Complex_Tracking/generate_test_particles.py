#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

# generate test particles for testing nll element.
# 2.5 MeV protons
# normalized emittance 0.3e-6 m-rad
# betax = betay = 0.65

mp = synergia.foundation.pconstants.mp
etot = 0.0025 + mp
normemit = 0.3e-6
betax = 0.65
betay = 0.65

npart = 100000

commxx = synergia.utils.Commxx()

refpart = synergia.foundation.Reference_particle(1, mp, etot)
bunch = synergia.bunch.Bunch(refpart, npart, 1.0e10, commxx)

emit = normemit * refpart.get_beta()*refpart.get_gamma()
sxsq = emit*betax
sxpsq = emit/betax
sysq = emit*betay
sypsq = emit/betay

means = np.zeros((6),'d')
corr = np.zeros((6,6),'d')
corr[0,0] = sxsq
corr[1,1] = sxpsq
corr[2,2] = sysq
corr[3,3] = sypsq

seed = 12345679
dist = synergia.foundation.Random_distribution(seed, commxx)
synergia.bunch.populate_transverse_gaussian(dist, bunch, means, corr, 0.01)

particles = bunch.get_local_particles()
np.savetxt("particles.txt", particles[:,0:4])

diag = synergia.bunch.Diagnostics_particles("particles.h5")
diag.set_bunch(bunch)
diag.update_and_write()
