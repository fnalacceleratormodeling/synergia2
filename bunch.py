#!/usr/bin/env python

import Numeric

from Pgrid2dPkgpy import *
from BeamBunchPkgpy import *
from ExternalBLEPkgpy import *
from DistributionPkgpy import *
from OutputPkgpy import *

from mpi4py import MPI

# Oi. We need a more general solution for this
MPI.MPI_COMM_WORLD=91

class Bunch:
    def __init__(self, current, beam_parameters, num_particles):
        self.beam_parameters = beam_parameters
        kinetic_energy = beam_parameters.kinetic_energy_GeV
        mass = beam_parameters.mass_GeV
        charge = beam_parameters.charge_e
        initial_phase = beam_parameters.initial_phase_rad
        self.beambunch = BeamBunch()
        construct_BeamBunch_external(self.beambunch, current,
                                     kinetic_energy * 1.0e9,
                                     mass* 1.0e9, charge, num_particles,
                                     initial_phase)
        nparam = 30
#         # obviously, I have punted here...
#         distparam = Numeric.array([ 0.000281124, 0, 4.9945e-07, 0, 0,
#                                     0.000991052, 0,
#                                     0, 0, 1.41675e-07, 0, 0, 0, 0,
#                                     0.345566, -6.62671e-07, 0, 0, 0, 0,
#                                     9.30938e-08, 0, 0, 0, 0, 0, 0,
#                                     1, 666, 666] ,'d')
        # punt number 2
        nprow = 1
        npcol = 1
        self.grid = Pgrid2d()
        construct_Pgrid2d_external(self.grid, MPI.MPI_COMM_WORLD, nprow, npcol)
        flagalloc = 0
        Gauss_Covariance_Dist_external(self.beambunch, nparam,
                                       beam_parameters.distparam(),
                                       self.grid, flagalloc)

    def current(self):
        return self.beambunch.Current
    def mass(self):
        return self.beambunch.Mass
    def charge(self):
        return self.beambunch.charge
    def num_particles(self):
        return self.beambunch.Npt
    def num_particles_local(self):
        return self.beambunch.Nptlocal
    def reference_particle(self):
        return self.beambunch.refptcl
    def particles(self):
        return self.beambunch.Pts1

    def write_fort(self,z):
        diagnostic3_Output(z, self.beambunch,
                           self.beam_parameters.scaling_frequency_Hz)
