#!/usr/bin/env python

import numpy
import synergia

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
    def __init__(self, ramp_turns, turns_to_extract, initial_k1, final_k1, final_k2l):
        synergia.simulation.Propagate_actions.__init__(self)
        Pickle_helper.__init__(self, ramp_turns, turns_to_extract, initial_k1, final_k1, final_k2l)
        self.ramp_turns = ramp_turns
        self.turns_to_extract = turns_to_extract
        self.initial_k1 = initial_k1
        self.final_k1 = final_k1
        self.final_k2l = final_k2l
    def turn_end_action(self, stepper, bunch, turn_num):
        synergia_elements = stepper.get_lattice_simulator().get_lattice().get_elements()
        turn_num += 1
        myrank = synergia.utils.Commxx().get_rank()
        # sextupole ramping
        if turn_num <= self.ramp_turns:
            index = 0
            length = 0
            for element in synergia_elements:
                if (turn_num == 50 or turn_num == 51):
                    length += element.get_length()
                    if myrank == 0:
                        print "%5d %10s %10g" % (turn_num, element.get_name(), length)
                if element.get_type() == "multipole":
                    new_k2l = self.final_k2l[index] * turn_num / self.ramp_turns
                    element.set_double_attribute("k2l", new_k2l)
                    index += 1
                    #if myrank == 0:
                    #    print
                    #    print "    turn                             :", 
                    #    print turn_num
                    #    print "    updated multipole                :", 
                    #    print element.get_name()
                    #    print "    final k2l                        :", 
                    #    print self.final_k2l[index - 1], "1/m^2"
                    #    print "    new k2l                          :", 
                    #    print new_k2l, "1/m^2"
                    #    print "    new_k2l (real)                   :", 
                    #    print element.get_double_attribute("k2l"), "1/m^2"
        # quadrupole ramping...
        if turn_num > self.ramp_turns:
            epsilon = 1.0 * (turn_num - self.ramp_turns) / self.turns_to_extract
            index = 0
            for element in synergia_elements:
                if element.get_type() == "quadrupole":
                    new_k1 = (1.0 - epsilon) * self.initial_k1[index] \
                            + epsilon * self.final_k1[index]
                    element.set_double_attribute("k1", new_k1)
                    index += 1
                    #if myrank == 0:
                    #    print
                    #    print "    turn                             :", 
                    #    print turn_num
                    #    print "    updated quadrupole               :", 
                    #    print element.get_name()
                    #    print "    epsilon                          :", epsilon
                    #    print "    initial k1                       :", 
                    #    print self.initial_k1[index - 1], "1/m"
                    #    print "    final k1                         :", 
                    #    print self.final_k1[index - 1], "1/m"
                    #    print "    new k1                           :", 
                    #    print new_k1, "1/m"
                    #    print "    new k1 (real)                    :", 
                    #    print element.get_double_attribute("k1"), "1/m"
        stepper.get_lattice_simulator().update()

