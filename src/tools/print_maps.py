#!/usr/bin/env python

import sys
import synergia
from mpi4py import MPI

class Options:
    def __init__(self):
        self.lattice_file = None
        self.lattice = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synprintmaps <filename> <lattice> [option1] ... [optionn]"
    print "available options are:"
    print "    no options yet!"
    sys.exit(0)

def print_maps(options):
    lattice = synergia.lattice.Mad8_reader().get_lattice(options.lattice,
                                                         options.lattice_file)
    for element in lattice.get_elements():
        element.set_string_attribute("extractor_type", "chef_map")
    map_order = 1
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice,
                                                              map_order)
    steps_per_element = 1
    stepper = synergia.simulation.Independent_stepper_elements(lattice_simulator,
                                                                steps_per_element)
    stepper.force_update_operations_no_collective()
    for step in stepper.get_steps():
        for operator in step.get_operators():
            if operator.get_type() == 'independent':
                io = synergia.simulation.as_independent_operator(operator)
                slices = io.get_slices()
                print "%s:" % slices[0].get_lattice_element().get_name()
                if len(slices) > 1:
                    raise RuntimeError, "found multiple slices in an independent operator"
                for operation in io.get_operations():
                    if operation.get_type() == 'fast_mapping':
                        fmo = synergia.simulation.as_fast_mapping_operation(operation)
                        dense_mapping = synergia.simulation.Dense_mapping(fmo.get_fast_mapping())
                        print "order 0:"
                        print dense_mapping.get_constant_term()
                        print "order 1:"
                        print dense_mapping.get_linear_term()
                        print

def handle_args(args):
    if len(args) < 2:
        do_help()
    options = Options()
    options.lattice_file = args[0]
    options.lattice = args[1]
    for arg in args[2:]:
        if arg == '--help':
            do_help(plotparams)
        else:
            do_error('Unknown argument "%s"' % arg)
    return options

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    print_maps(options)
