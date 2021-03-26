#!/usr/bin/env python

import synergia
import numpy as np

class Pickle_helper:
    __getstate_manages_dict__ = 1
    def __init__(self, *args):
        self.args = args
    def __getinitargs__(self):
        return self.args
    def __getstate__(self):
        return self.__dict__
    def __setstate__(self, state):
        self.__dict__ = state

class Ramp_actions(synergia.simulation.Propagate_actions, Pickle_helper):
    def __init__(self):
        synergia.simulation.Propagate_actions.__init__(self)
        fh = open("nf.dat", "w")
        self.fh = fh
        Pickle_helper.__init__(self)
    def turn_end_action(self, stepper, bunch, turn_num):
        particles = np.array(bunch.get_local_particles())
        nx, ny = particles.shape
        tmppart = np.zeros((nx,ny), dtype='d')
        tmppart[:,:] = particles[:,:]

        #print('Turn ', turn_num, ' num part: ', bunch.get_local_num())
        lattice_simulator = stepper.get_lattice_simulator()
        lattice_simulator.convert_xyz_to_normal(tmppart)
        partout = np.zeros(80 * 6)
        for p in range(bunch.get_local_num()):
            # some particles might have been lost so build up particles
            # to write out based on particle id
            id = int(tmppart[p, 6])
            for j in range(6):
                partout[id*6+j] = tmppart[p, j]
        for j in range(partout.shape[0]):
            if j != 0:
                print(" ", file=self.fh, end="")
            print(partout[j], file=self.fh, end="")
        print(file=self.fh)
        self.fh.flush()
