#!/usr/bin/env python

import sys
import time
import numpy
import nlopt
import synergia
from mpi4py import MPI
from plasma_options import opts

def update_triplet(simulator, x):
    for element in simulator.get_lattice().get_elements():
        name = element.get_name()
        type = element.get_type()
        if name == "modq4d":
            element.set_double_attribute("k1", x[0])
        elif name == "modq3f":
            element.set_double_attribute("k1", x[1])
        elif name == "modq2d":
            element.set_double_attribute("k1", x[2])
    simulator.update()

def tracking(x):
    #lattice_Tev = synergia.lattice.Mad8_reader().get_lattice("injection", "tevatron.lat")
    #lattice_simulator_Tev = synergia.simulation.Lattice_simulator(lattice_Tev, opts.map_order)
    global lattice_simulator_Tev
    #global bunch
    bunch = synergia.optics.generate_matched_bunch_transverse(
              lattice_simulator_Tev, opts.emit, opts.emit, opts.stdz,
              opts.dpop, opts.real_particles, opts.macro_particles,
              seed=opts.seed)

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

    #lattice_line = synergia.lattice.Mad8_reader().get_lattice("plasma", "tevatron.lat")
    #lattice_simulator = synergia.simulation.Lattice_simulator(lattice_line, opts.map_order)

    global lattice_simulator
    update_triplet(lattice_simulator, x)

    global count
    stepper = synergia.simulation.Independent_stepper_elements(
                                lattice_simulator, opts.steps)
    diagnostics_full2 = synergia.bunch.Diagnostics_full2("full2_%02d.h5" % count)   
    diagnostics_particles = synergia.bunch.Diagnostics_particles("particles_%02d.h5" % count)
    bunch_simulator.add_per_step(diagnostics_full2)
    bunch_simulator.add_per_turn(diagnostics_particles)

    propagator = synergia.simulation.Propagator(stepper)
    propagator.propagate(bunch_simulator, opts.turns)

    x_std = diagnostics_full2.get_std()[0]
    y_std = diagnostics_full2.get_std()[2]

    #diagnostics_step = synergia.bunch.Diagnostics_full2(bunch, ("full2_%03d.h5" % count))
    #diagnostics_turn = synergia.bunch.Diagnostics_particles(bunch, ("particles_%03d.h5" % count))

    #propagator = synergia.simulation.Propagator(stepper)
    #propagator.propagate(bunch, opts.turns, diagnostics_step, diagnostics_turn, opts.verbose)

    #x_std = diagnostics_step.get_std()[0]
    #y_std = diagnostics_step.get_std()[2]

    return x_std, y_std

def twiss(x):
    lattice_sextant=synergia.lattice.Mad8_reader().get_lattice("lo_beta_sextant", "tevatron.lat")
    sextant_elements = lattice_sextant.get_elements()
    sextant_simulator = synergia.simulation.Lattice_simulator(lattice_sextant, opts.map_order)
    update_triplet(sextant_simulator, x)
    global count
    if rank == 0:
        twiss_file = ("twiss_%03d.txt" % count)
        twiss_log = open(twiss_file, "w")
    index = 0
    length = 0.0
    for element in sextant_elements:
        lattice_functions = sextant_simulator.get_lattice_functions(element)
        type = element.get_type()
        name = element.get_name()
        length += element.get_length()
        beta_x = lattice_functions.beta_x
        beta_y = lattice_functions.beta_y
        alpha_x = lattice_functions.alpha_x
        D_x = lattice_functions.D_x
        if rank == 0:
            twiss_log.write("%3d %12s %10s %10.5f %10.5f %10.5f %10.5f %10.5f\n" % (
                        index, type, name, length, beta_x, beta_y, alpha_x, D_x))
            twiss_log.flush()
        index += 1
    if rank == 0:
        twiss_log.close()

def myfunc(x, grad):
    global count
    retval = numpy.zeros([6], 'd')
    retval = tracking(x)
    count += 1
    if rank == 0:
        print "Try %2g: (%10.7g, %10.7g, %10.7g) = %10.7g, %10.7g" % (count, \
                x[0], x[1], x[2], retval[0], retval[1])

    return retval[0] * retval[0] + retval[1] * retval[1]

def myconstraint1(x, grad):
    #retval = tracking(x)[0] - 0.0050
    retval = tracking(x)[0] - 0.00557097
    return retval

def myconstraint2(x, grad):
    #retval = tracking(x)[1] - 0.0050
    retval = tracking(x)[1] - 0.0049789
    return retval

###############################################################################

t0 = time.time()
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

#   Initial bunch and lattice_simulator setups
lattice_Tev = synergia.lattice.Mad8_reader().get_lattice("injection", "tevatron.lat")
lattice_simulator_Tev = synergia.simulation.Lattice_simulator(lattice_Tev, opts.map_order)

#~bunch = synergia.optics.generate_matched_bunch_transverse(
#~                  lattice_simulator_Tev, opts.emit, opts.emit, opts.stdz,
#~                  opts.dpop, opts.real_particles, opts.macro_particles,
#~                  seed=opts.seed)

lattice_line = synergia.lattice.Mad8_reader().get_lattice("plasma", "tevatron.lat")
lattice_simulator = synergia.simulation.Lattice_simulator(lattice_line, opts.map_order)

count = 0
#algorithm = nlopt.LD_MMA
algorithm = nlopt.LN_COBYLA
opt = nlopt.opt(algorithm, 3)
if rank == 0:
    print
    print "Algorithm:", opt.get_algorithm_name()

#   Ojbective function
opt.set_min_objective(myfunc)

#   Bound constraints
lb = numpy.array([-0.06, 0.01, -0.06], 'd')
ub = numpy.array([-0.01, 0.06, -0.01], 'd')
opt.set_lower_bounds(lb)
opt.set_upper_bounds(ub)

#   Nonlinear constraints
opt.add_inequality_constraint(myconstraint1, 1.0e-4)
opt.add_inequality_constraint(myconstraint2, 1.0e-4)

#   Stopping criteria
tol = 1.0e-4
#Set relative on optimization parameters.
opt.set_xtol_rel(tol)
#opt.set_ftol_rel(tol)

#   Performing the optimization
#minuit's answer
x = numpy.array([-2.898641E-02, 3.712463E-02, -4.090105E-02], 'd')
#x = numpy.array([-0.040060272135406534, 0.040026874437379276, -0.040060272135406534], 'd')
xopt = opt.optimize(x)

opt_val = opt.last_optimum_value()
result = opt.last_optimize_result()

if rank == 0:
    print
    print "found optimum at x[0] = %10.8f, x[1] = %10.8f, and x[2] = %10.8f" % (xopt[0], xopt[1], xopt[2])
    print "after %03d evaluations, x_std = %10.8f" % (count, opt_val)

t1 = time.time()

if rank == 0:
    print "elapsed time =", t1 - t0
