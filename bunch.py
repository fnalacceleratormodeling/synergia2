#!/usr/bin/env python

import Numeric

from Pgrid2dPkgpy import *
from BeamBunchPkgpy import *
from ExternalBLEPkgpy import *
from DistributionPkgpy import *
from OutputPkgpy import *

import loadfile

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
    def generate_particles(self):
        flagalloc = 0
        print "transvers =",self.beam_parameters.get_transverse()
        if self.beam_parameters.get_transverse():
            print "generating transverse"
            Gauss_Covariance_Trans_Dist_external(\
                self.beambunch,
                self.beam_parameters.get_nparam(),
                self.beam_parameters.get_distparam(),
                self.processor_grid.get_pgrid2d(),
                flagalloc)
        else:
            print "generating not transverse"
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
    
    def write_fort(self,z):
        diagnostic3_Output(z, self.beambunch,
                           self.beam_parameters.scaling_frequency_Hz)
    def write_particles(self,filename):
        f = open(filename,"w")
        p = self.particles()
        for i in range(0,p.shape[1]):
#            f.write("%g %g\n" % tuple([p[0,i],p[1,i]))
            f.write("%g %g %g %g %g %g %g\n" % tuple(p[:,i].tolist()))
        f.close()
        
