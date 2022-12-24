#!/usr/bin/env python
import sys
import os
import numpy as np
import synergia
import synergia.simulation as SIM
ET = synergia.lattice.element_type
MT = synergia.lattice.marker_type
#####################################


#######################################################



#######################################################


################################################################################



################################################################################


################################################################################


################################################################################

def get_lattice():
    # read the lattice in from a MadX sequence file
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "cfoborodobo32.madx")
    lattice.set_all_string_attribute("extractor_type", "libff")
    return lattice

################################################################################

def make_stepper():
    # space charge
    if opts.spacecharge:

        if opts.spacecharge == "3dopen-hockney":
            sc_ops = synergia.collective.Space_charge_3d_open_hockney_options(gridx, gridy, gridz)
            sc_ops.comm_group_size = opts.comm_group_size

        elif solver == "2dopen-hockney":
            sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(gridx, gridy, gridz)
            sc_ops.comm_group_size = opts.comm_group_size

        elif solver == "rectangular":
            space_charge = synergia.collective.Space_charge_rectangular(grid, opts.pipe_size)
            sc_ops.comm_group_size = opts.comm_group_size

        else:
            sys.stderr.write("foborodobo32.py: solver must be either 3dopen-hockney, 2dopen-hockney, or rectangular\n")
            sys.exit(1)

        stepper = synergia.simulation.Split_operator_stepper(sc_ops, opts.steps)

    else:
        # space charge not active

        if opts.stepper == "splitoperator":
            sc_ops = synergia.collective.Dummy_CO_options()
            stepper = synergia.simulation.Split_operator_stepper(sc_ops, opts.steps)

        elif opts.stepper == "elements":
            stepper = synergia.simulation.Independent_stepper_elements(opts.steps)

        else:
            sys.stderr.write("mi.py: stepper must be either splitopertor,independent, or elements\n")
            sys.exit(1)
    return stepper


################################################################################

class opts:
    pass

def init():
    opts.energy = 0.8
    opts.comm_group_size = 16
    opts.spacecharge = False
    #opts.stepper = "splitoperator"
    #opts.steps = 96
    opts.stepper = "elements"
    opts.steps = 1


################################################################################

def main():

    init()

    logger = synergia.utils.Logger(0)
    lattice = get_lattice()
    print('Read lattice, length = {}, {} elements'.format(lattice.get_length(), len(lattice.get_elements())), file=logger)

    # set the momentum of the beam.  This could have been in a beam statement
    # in the lattice file, but this gives the option of setting it on the
    # command line by creating a reference particle with the desired momentum.

    # create the Reference particle object with the correct energy
    total_energy = opts.energy + synergia.foundation.pconstants.mp

    refpart = synergia.foundation.Reference_particle(1, synergia.foundation.pconstants.mp, total_energy)

    # set it into the lattice object
    lattice.set_reference_particle(refpart)
    energy = refpart.get_total_energy()
    momentum = refpart.get_momentum()
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    print("energy: ", energy, file=logger)
    print("momentum: ", momentum, file=logger)
    print("gamma: ", gamma, file=logger)
    print("beta: ", beta, file=logger)

    orbit_length = lattice.get_length()
    print('Orbit length: ', orbit_length, file=logger)

    f = open("cfoborodobo32_lattice.out", "w")
    print(lattice, file=f)
    f.close()

    stepper = make_stepper()

    # propagator
    propagator = synergia.simulation.Propagator(lattice, stepper)

    # logs
    #simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_STEP)
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN)
    screen = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.DEBUG)

    print('number of lattice elements: ', len(lattice.get_elements()))
    print('lattice first element a1 param: ', lattice.get_elements()[0].get_double_attribute('a1'))
    print()
    lattice.get_elements()[0].set_double_attribute('a1', -9999.0)
    print('lattice first element a1 param after mod: ', lattice.get_elements()[0].get_double_attribute('a1'))
    print('Modifying the lattice element directly from the lattice works')
    print()
    print()
    print('len(propagator.get_lattice().get_elements()): ', len(propagator.get_lattice().get_elements()))


    proplattice = propagator.get_lattice()
    print('number propagator.get_lattice(): ', len(proplattice.get_elements()))
    print('number of propagator.get_lattice_elements: ', len(propagator.get_lattice_elements()))
    print()
    # print()
    # print('dir(lattice_mutable): ', dir(lattice_mutable))
    # print()
    # print('dir(lattice_const): ', dir(lattice_const))
    print()
    print('propagator.get_lattice_elements first element a1 param orig: ', propagator.get_lattice_elements()[0].get_double_attribute('a1'))
    print()
    propagator.get_lattice_elements()[0].set_double_attribute('a1', 9999.0)
    print('propagator.get_lattice_elements first element a1 after nod mod: ', propagator.get_lattice_elements()[0].get_double_attribute('a1'))
    print('Modifying the lattice element from propagator.get_lattice does NOT work')
    print()
    print()
    print('proplattice first element a1 param orig: ', proplattice.get_elements()[0].get_double_attribute('a1'))
    proplattice.get_elements()[0].set_double_attribute('a1', -222.0)
    print('proplattice first element mod: ', proplattice.get_elements()[0].get_double_attribute('a1'))
    print('Modifying the lattice element from propagator.get_lattice().get_elements() appears to work, or at tleast the readback shows the change')
    i=0
    slices = propagator.get_lattice_element_slices()
    for s in slices:
        if i >= 1:
            break
        print(s)
        i = i+1
    print(dir(list(slices)[0]))
    #list(slices)[0].print_()
    print(list(slices)[0].get_lattice_element())
    print('The change to a1 is reflected in the lattice slices.')

    proplattice.get_reference_particle().set_total_energy(100.0)
    print('modified reference particle energy: ', proplattice.get_reference_particle().get_total_energy())

if __name__ == "__main__":
    main()

