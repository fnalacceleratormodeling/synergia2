#!/usr/bin/env python

import sys
import numpy
import nlopt
from numpy import *

def myfunc(x, grad):
    global count
    count += 1
    if grad.size > 0:
        #grad is a numpy array of length n, which should be set to then 
        #gradient of the function wi.r.t the optimization parameters at x.
        grad[0] = 0.0
        grad[1] = 0.5 / numpy.sqrt(x[1])
    print "try %2d: x[0] = %10.8f, x[1] = %10.8f, f(x[0], x[1]) = %10.8f" % (count, x[0], x[1], sqrt(x[1]))
    return numpy.sqrt(x[1])

def get_grad_constraint(x, data):
    grad = numpy.array([0,0], 'd')
    grad[0] = 3 * data[0] * (data[0] * x[0]) * (data[0] * x[0] + data[1])
    grad[1] = -1.0
    return grad

def get_constraint(x, data):
    return ((data[0] * x[0] + data[1]) * (data[0] * x[0] + data[1]) \
            * (data[0] * x[0] + data[1]) - x[1])

def myconstraint1(x, grad):
    data = numpy.array([2,0], 'd')
    retval = get_constraint(x, data)
    if grad.size > 0:
        grad = get_grad_constraint(x, data)
    return retval

def myconstraint2(x, grad):
    data = numpy.array([-1,1], 'd')
    retval = get_constraint(x, data)
    if grad.size > 0:
        grad = get_grad_constraint(x, data)
    return retval


#   The nlopt.opt class

#algorithm = nlopt.LD_MMA        # Method of Moving Asymptotes (local, derivative)
algorithm = nlopt.LN_COBYLA    # Constrained Optimization BY Linear Approximations (local, non-derivative)
n = 2                           # the number of optimization parameters
opt = nlopt.opt(algorithm, 2)
print opt.get_algorithm_name(), opt.get_dimension()


#   Ojbective function

opt.set_min_objective(myfunc)
#opt.set_max_objective(myfunc)


#   Bound constraints

lb = numpy.array([0,0], 'd')
lb[0] = - numpy.inf             # an unbounded dimension
lb[1] = 0.0
opt.set_lower_bounds(lb)
#print opt.get_lower_bounds()

#opt.set_upper_bounds(ub)
#print opt.get_upper_bounds()


#   Nonlinear constraints

opt.add_inequality_constraint(myconstraint1, 1.0e-8)
opt.add_inequality_constraint(myconstraint2, 1.0e-8)
#opt.remove_inequality_constraints()

#opt.add_equality_constraint(myconstraint, tol=1.0e-8)
#opt.remove_equality_constraint()


#   Stopping criteria

tol = 1.0e-8
#Stop when an objective value of at least stopval is found.
#opt.set_stopval(stopval)
#opt.get-stopval()

#Set relative on function value.
opt.set_ftol_rel(tol)
opt.get_ftol_rel()

#Set absolute tolerance on function value.
#opt.set_ftol_abs(tol)
#opt.get_ftol_abs()

#Set relative on optimization parameters.
#opt.set_xtol_rel(tol)
#opt.get_xtol_rel()

#Set absolute tolerance on optimization parameters.
#opt.set_xtol_abs(tol)
#opt.get_xtol_abs()

#Stop when the number of funtion evaluations exceeds maxeval.
#opt.set_maxeval(maxeval)
#opt.get_maxeval()

#Stop when the optimization time (in secs) excees maxtime.
#opt.set_maxtime(maxtime)
#opt.get_maxtime()

#nlopt.Forcesstop


#   Performing the optimization

count = 0
x = numpy.array([1.234, 5.678], 'd')
xopt = opt.optimize(x)

opt_val = opt.last_optimum_value()
result = opt.last_optimize_result()

print "found optimum at x[0] = %10.8f, x[1] = %10.8f" % (xopt[0], xopt[1])
print "after %2d evaluations, f(x[0], x[1]) = %10.8f" % (count, opt_val)

