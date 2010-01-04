#!/usr/bin/env python

import synergia

from bmlfactory import *
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *
from synergia.physics_constants import *

import debuncher

import math
import sys
import numpy
import string
import os.path

class Dgourmet(synergia.Gourmet):
    def __init__(self, initial_kinetic_energy,
                 scaling_frequency, order,tune_h,tune_v,write_output=True, particle='proton'):
        self.write_output = write_output
        if self.write_output:
            filename = "chef_log.txt"
        else:
            filename = "/dev/null"
        self.debuncher = debuncher.Debuncher(filename,tune_h,tune_v)
        #~ self.lattice_file = lattice_file
        #~ self.line_name = line_name
        self.scaling_frequency = scaling_frequency
        self.order = order
        self.saved_elements = []
        JetParticle.createStandardEnvironments(self.order)

        self.particle = particle
        if self.particle == 'proton':
            self.mass = PH_NORM_mp
        elif self.particle == 'positron':
            self.mass = PH_NORM_me
        else:
            raise RuntimeError, 'Unknown particle %s' % self.particle
        self.initial_kinetic_energy = initial_kinetic_energy
        self.initial_energy = self.initial_kinetic_energy + self.mass
        self.final_energy = None
        self.initial_momentum = math.sqrt(self.initial_energy**2 -
                                          self.mass**2)
        
        self.beamline = self.debuncher.get_beamline()
        
        
    def complete_setup(self):
        self.debuncher.complete_setup()
        JetParticle.createStandardEnvironments(self.order)
        self.beamline = self.debuncher.get_beamline()
        self.needs_commission = False
        self.is_commissioned = False
        ##self._commission()
        self.have_actions = 0
        self.have_fast_mappings = 0
        self.have_element_fast_mappings = 0
        self.context = BeamlineContext(self.get_initial_particle(),
                                       self.beamline)
        if not self.context.isTreatedAsRing():
            self.context.handleAsRing()

    def get_sextupoles(self):
        values = []
        for element in self.beamline:
            if element.Type() == 'thinSextupole':
                values.append(element.Strength())
        return numpy.array(values,'d')
        
    def set_sextupoles(self,values):
        index = 0
        for element in self.beamline:
            if element.Type() == 'thinSextupole':
                print "setting type of",element.Name(),"to",values[index]
                element.setStrength(values[index])
                index += 1
                
        