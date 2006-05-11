#!/usr/bin/env python

import Numeric

from UberPkgpy import *
from OutputPkgpy import *
from mpi4py import MPI

import loadfile

# this is a comment

seed_offset = 0
class Bunch:
    def __init__(self, current, beam_parameters, num_particles,
                 processor_grid):
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
        self.processor_grid = processor_grid
        global seed_offset
        init_seed_external(processor_grid.get_pgrid2d(),seed_offset)
        seed_offset += MPI.size
        
    def generate_particles(self):
        flagalloc = 0
        if self.beam_parameters.get_transverse():
            Gauss_Covariance_Trans_Dist_external(\
                self.beambunch,
                self.beam_parameters.get_nparam(),
                self.beam_parameters.get_distparam(),
                self.processor_grid.get_pgrid2d(),
                flagalloc)
        else:
            Gauss_Covariance_Dist_external(\
                self.beambunch,
                self.beam_parameters.get_nparam(),
                self.beam_parameters.get_distparam(),
                self.processor_grid.get_pgrid2d(),
                flagalloc)

    def read_particles(self,filename):
        self.beambunch.Pts1 = loadfile.loadfile_transpose(filename)
        self.beambunch.Npt = self.beambunch.Pts1.shape[1]
        self.beambunch.Nptlocal = self.beambunch.Npt
        print "read",self.beambunch.Npt,"particles"
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
    def get_beambunch(self):
        return self.beambunch
    def inject(self,injected):
        inject_BeamBunch_external(self.get_beambunch(),
                                  injected.get_beambunch())
    def write_fort(self,z):
        diagnostic3_Output(z, self.beambunch,
                           self.beam_parameters.scaling_frequency_Hz)
    def write_particles(self,filename):
        if MPI.rank == 0:
            f = open(filename,"w")
            for proc in xrange(1,MPI.size):
                parts = MPI.WORLD.Recv(source=proc)
                for i in range(0,parts.shape[1]):
                    f.write("%g %g %g %g %g %g %g\n" % \
                            tuple(parts[:,i]))
            parts = self.particles()
            for i in range(0,parts.shape[1]):
                f.write("%g %g %g %g %g %g %g\n" % \
                        tuple(parts[:,i]))
            f.close()
        else:
            MPI.WORLD.Send(self.particles(),dest=0)
