#!/usr/bin/env python3

import sys, os
from mpi4py import MPI
import numpy as np
import synergia

from aperture_demo_options import opts

macroparticles = opts.macroparticles
real_particles = opts.real_particles
turns = 1

def get_lattice():

    fodo_madx = """
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1= 0.071428571428571425;
d: quadrupole, l=2.0, k1= -0.071428571428571425;
m: marker;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
ap1: m, at = 2.0;
fodo_2: o, at=2.0;
ap2: m at = 10.0;
fodo_3: d, at=10.0;
ap3: m at = 12.0;
fodo_4: o, at=12.0;
ap4 : m at = 20.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice('fodo')

    return lattice

def print_statistics(bunch):

    parts = bunch.get_particles_numpy()
    print(parts.shape,  ", ", parts.size )
    print("shape: {0}, {1}".format(parts.shape[0], parts.shape[1]))

    mean = synergia.bunch.Core_diagnostics.calculate_mean(bunch)
    std = synergia.bunch.Core_diagnostics.calculate_std(bunch, mean)
    print("mean = {}".format(mean))
    print("std = {}".format(std))

def create_simulator(ref_part):

    print('in create_simulator')
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, macroparticles, real_particles)
    print('simlator created with single bunch')

    bunch = sim.get_bunch()

    bunch_means = np.zeros(6, dtype='d')
    # numerical values from cxx_covariance_matrix.xml in cxx_benchmark
    bunch_covariances = np.array(
        [[3.0509743977035345e-05, 2.2014134466660509e-06, 0, 0, 0, 0],
         [2.2014134466660509e-06, 1.9161816525115869e-07, 0, 0, 0, 0],
         [0, 0, 7.5506914064526925e-06, -6.6846812465678249e-07, 0, 0],
         [0, 0, -6.6846812465678249e-07, 1.9161816525115867e-07, 0, 0],
         [0, 0, 0, 0, 0.00016427607645871527, 0],
         [0, 0, 0, 0, 0, 1e-08]])
    #print('bunch_covariances.shape: ', bunch_covariances.shape)
    
    dist = synergia.foundation.PCG_random_distribution(1234567, synergia.utils.Commxx())

    synergia.bunch.populate_6d( dist, 
        bunch, 
        bunch_means,
        bunch_covariances)


    # zero 0 momenta to make particles go straight forward
    parts = bunch.get_particles_numpy()
    parts[:,1] = 0.0
    parts[:,3] = 0.0
    

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Dummy_CO_options()

    stepper = synergia.simulation.Split_operator_stepper_elements(sc_ops, 1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

class mydiag(synergia.bunch.Diagnostics):
    def __init__(self, filename):
        synergia.bunch.Diagnostics.__init__(self, "mydiag", filename)
        print("mydiag created")

    def do_update(self, bunch):
        print("my diag update")

    def do_reduce(self, comm, root):
        print("my diag reduce")

def run():

    screen = synergia.utils.parallel_utils.Logger(0, 
            #synergia.utils.parallel_utils.LoggerV.DEBUG)
                                                  synergia.utils.parallel_utils.LoggerV.INFO)

    lattice = get_lattice()
    print('Read lattice, length: ', lattice.get_length(), ', ', len(lattice.get_elements()), ' elements', file=screen)
    

    if len(lattice.get_elements()) != 8:
        print('ERROR!!! Wrong number of lattice elements!!!!!')
        for elem in lattice.get_elements():
            elem.print_()

    # find the aperture screen
    ap = None
    for elem in lattice.get_elements():
        if elem.get_name()[0:2] == "ap":

            if opts.aperture:
                print('setting aperture type ', opts.aperture, file=screen)
            if opts.hoffset:
                print('    horizontal offset: ', opts.hoffset, file=screen)
            if opts.voffset:
                print('    vertical offset: ', opts.voffset, file=screen)

            if opts.aperture == "circular":
                print('circular aperture radius: ', opts.circular_aperture_radius)
                elem.set_string_attribute("aperture_type","circular")
                elem.set_double_attribute("circular_aperture_radius", opts.circular_aperture_radius)
                if opts.hoffset:
                    elem.set_double_attribute('hoffset', opts.hoffset)
                if opts.voffset:
                    elem.set_double_attribute('voffset', opts.voffset)
            elif opts.aperture == "elliptical":
                elem.set_string_attribute("aperture_type","elliptical")
                elem.set_double_attribute("elliptical_aperture_horizontal_radius", opts.elliptical_horizontal_radius)
                elem.set_double_attribute("elliptical_aperture_vertical_radius", opts.elliptical_vertical_radius)
                if opts.hoffset:
                    elem.set_double_attribute('hoffset', opts.hoffset)
                if opts.voffset:
                    elem.set_double_attribute('voffset', opts.voffset)
            elif opts.aperture == "rectangular":
                elem.set_string_attribute("aperture_type","rectangular")
                elem.set_double_attribute("rectangular_aperture_width", 2*0.005)
                elem.set_double_attribute("rectangular_aperture_height", 2*0.002)
                if opts.hoffset:
                    elem.set_double_attribute('hoffset', opts.hoffset)
                if opts.voffset:
                    elem.set_double_attribute('voffset', opts.voffset)
            elif opts.aperture == "polygon":
                elem.set_string_attribute("aperture_type","polygon")
                elem.set_double_attribute("the_number_of_vertices", 4)

                elem.set_double_attribute("pax1", opts.polygon_aperture_arm)
                elem.set_double_attribute("pay1", 0.0)
                elem.set_double_attribute("pax2", 0.0)
                elem.set_double_attribute("pay2", opts.polygon_aperture_arm)
                elem.set_double_attribute("pax3", -opts.polygon_aperture_arm)
                elem.set_double_attribute("pay3", 0.0)
                elem.set_double_attribute("pax4", 0.0)
                elem.set_double_attribute("pay4", -opts.polygon_aperture_arm)
                elem.set_double_attribute('hoffset', opts.hoffset)
                elem.set_double_attribute('voffset', opts.voffset)
            else:
                print("unknown aperture type '%s'" % opts.aperture)

        if elem.get_type_name() == "quadrupole":
            elem.set_double_attribute('k1', 0.0)

    sim = create_simulator(lattice.get_reference_particle())

    propagator = create_propagator(lattice)

    for s in propagator.get_lattice_element_slices():
        print(s)

    class context:
        steps = 0

    def action(sim, lattice, turn, step):
        #nonlocal another_steps
        #another_steps += 1
        context.steps += 1

    sim.reg_prop_action_step_end(action)

    # diagnostics

    diag_part = synergia.bunch.Diagnostics_particles("particles.h5")
    sim.reg_diag_per_step(diag_part)

    # logger
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_TURN)
            #synergia.utils.parallel_utils.LoggerV.INFO)
            #synergia.utils.parallel_utils.LoggerV.INFO_STEP)

    # propagate
    propagator.propagate(sim, simlog, turns)

    print("total steps = ", context.steps)

    print('after propagate')
    sim.get_bunch().print_statistics(screen);

    # save
    #synergia.simulation.checkpoint_save(propagator, sim);

def main():

    try:
        run()
    except:
        raise RuntimeError("Failure to launch aperture_demo")

main()
