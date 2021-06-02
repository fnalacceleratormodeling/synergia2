#!/usr/bin/env python3

import sys, os
from mpi4py import MPI
import numpy as np
import synergia

macroparticles = 524288
real_particles=2.94e10

def get_lattice():

    fodo_madx = """
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1=0.071428571428571425;
d: quadrupole, l=2.0, k1=-0.071428571428571425;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
fodo_2: o, at=2.0;
fodo_3: d, at=10.0;
fodo_4: o, at=12.0;
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
    #dist = synergia.foundation.Random_distribution(1234567, synergia.utils.Commxx())

    #print('before populate_6d')
    #print_statistics(bunch)
    synergia.bunch.populate_6d( dist, 
        bunch, 
        bunch_means,
        bunch_covariances)

    #print('after populate_6d')
    #print_statistics(bunch)

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Space_charge_2d_open_hockney_options(64, 64, 64)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 4)
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
            synergia.utils.parallel_utils.LoggerV.DEBUG)

    lattice = get_lattice()
    print('Read lattice, length: ', lattice.get_length(), ', ', len(lattice.get_elements()), ' elements', file=screen)
    if len(lattice.get_elements()) != 4:
        print('ERROR!!! Wrong number of lattice elements!!!!!')
        for elem in lattice.get_elements():
            elem.print_()
    sim = create_simulator(lattice.get_reference_particle())

    propagator = create_propagator(lattice)

    sim.get_bunch().print_statistics(screen);

    class context:
        steps = 0

    def action(sim, lattice, turn, step):
        #nonlocal another_steps
        #another_steps += 1
        context.steps += 1

    sim.reg_prop_action_step_end(action)

    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    sim.reg_diag_per_turn(diag_full2)

    diag_bt = synergia.bunch.Diagnostics_bulk_track("diag_track.h5", 1000, 0)
    sim.reg_diag_per_turn(diag_bt)

    diag_part = synergia.bunch.Diagnostics_particles("diag_part.h5", 100)
    sim.reg_diag_per_turn(diag_part)

    sim.reg_diag_per_turn(mydiag("mydiag.h5"))
    diag_dummy = synergia.bunch.Diagnostics_dummy()
    sim.reg_diag_per_turn(diag_dummy)

    # logger
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_STEP)

    # propagate
    propagator.propagate(sim, simlog, 10)

    print("total steps = ", context.steps)

    print('after propagate')
    sim.get_bunch().print_statistics(screen);

    # save
    #synergia.simulation.checkpoint_save(propagator, sim);

def main():

    print("running fodo.py: my rank =", MPI.COMM_WORLD.Get_rank())
    try:
        run()
    except:
        raise RuntimeError("Failure to launch fodo.run")

main()
