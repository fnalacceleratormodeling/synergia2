#!/usr/bin/env python

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
    # The arguments to __init__ are what the Ramp_actions instance is
    # initialized with
    def __init__(self, multiplier):
        synergia.simulation.Propagate_actions.__init__(self)
        # pickling the arguments to the initializer allows the
        # module to resume after checkpointing.  They should be in the same
        # order as the arguments to __init__.
        Pickle_helper.__init__(self, multiplier)
        self.multiplier = multiplier
    def turn_end_action(self, stepper, bunch, turn_num):
        print "modifying lattice"
        for element in stepper.get_lattice_simulator().get_lattice().get_elements():
            if element.get_type() == "quadrupole":
                old_k1 = element.get_double_attribute("k1")
                element.set_double_attribute("k1", self.multiplier*old_k1)
                print "  updated", element.get_name(),"k1 =", self.multiplier*old_k1
        stepper.get_lattice_simulator().update()


