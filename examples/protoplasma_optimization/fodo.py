#!/usr/bin/env python

import sys
import math
import numpy
import nlopt
import synergia
from mpi4py import MPI
from fodo_options import opts

def update_fodo(simulator, x):
    for element in simulator.get_lattice().get_elements():
        name = element.get_name()
        type = element.get_type()
        if name == "f":
            element.set_double_attribute("k1", x[0])
        elif name == "d":
            element.set_double_attribute("k1", x[1])
    simulator.update()

def fodo(x):
    lattice = synergia.lattice.Mad8_reader().get_lattice("fodo", "fodo.lat")
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
    bunch = synergia.optics.generate_matched_bunch_transverse(
                  lattice_simulator, opts.emit, opts.emit, opts.stdz, opts.dpop,
                  opts.real_particles, opts.macro_particles,
                  seed=opts.seed)

    bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
    update_fodo(lattice_simulator, x)

    global count
    stepper = synergia.simulation.Independent_stepper_elements(
                                lattice_simulator, opts.steps)


    diagnostics_full2 = synergia.bunch.Diagnostics_full2("full2_%02d.h5" % count)
    diagnostics_particles = synergia.bunch.Diagnostics_particles("particles_%02d.h5" % count)
    bunch_simulator.add_per_step(diagnostics_full2)
    bunch_simulator.add_per_turn(diagnostics_particles)

    propagator = synergia.simulation.Propagator(stepper)
    #propagator.set_checkpoint_period(opts.checkpointperiod)
    propagator.propagate(bunch_simulator, opts.turns)

    x_std = diagnostics_full2.get_std()[0]
    y_std = diagnostics_full2.get_std()[2]
    return x_std, y_std

def myfunc(x, grad):
    global count
    retval = numpy.zeros([6], 'd')
    retval = fodo(x)
    spotsize = (retval[0]-0.010)**2
    count += 1
    if rank == 0:
        print "Try %2g: (%10.7g, %10.7g) = %10.7g, %10.7g, %10.7g" % (count, \
                x[0], x[1], retval[0], retval[1], spotsize)

    return spotsize

def myconstraint(x, grad):
# Funky way to constrain both strengths to be less than 0.5
    retval = 1.0
    if (math.fabs(x[0])-0.5) < 0.0:
        if (math.fabs(x[1])-0.5)< 0.0:
            retval=-1.0
    return retval

###############################################################################

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
count = 0
algorithm = nlopt.LN_COBYLA
opt = nlopt.opt(algorithm, 2)

if rank == 0:
    print
    print "Algorithm:", opt.get_algorithm_name()

#   Ojbective function
opt.set_min_objective(myfunc)

#   Bound constraints
lb = numpy.array([0.01, -0.15], 'd')
ub = numpy.array([0.15, -0.01], 'd')
opt.set_lower_bounds(lb)
opt.set_upper_bounds(ub)

#   Nonlinear constraints
opt.add_inequality_constraint(myconstraint, 1.0e-4)

#   Stopping criteria
tol = 1.0e-4
opt.set_xtol_rel(tol)
#opt.set_ftol_rel(tol)

#   Performing the optimization
focus = 7.0
length = 2.0
strength = 1.0 /(focus * length)
x = numpy.array([strength, -strength], 'd')
xopt = opt.optimize(x)

opt_val = opt.last_optimum_value()
result = opt.last_optimize_result()

if rank == 0:
    print
    print "found optimum at x[0] = %10.8f, x[1] = %10.8f" % (xopt[0], xopt[1])
    print "after %02d evaluations, (x_std-0.01)**2 = %10.8f" % (count, opt_val)
