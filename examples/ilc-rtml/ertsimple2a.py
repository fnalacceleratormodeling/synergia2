#!/usr/bin/env python

import local_paths

from bmlfactory import *
from basic_toolkit import *
from mxyzptlk import *
from beamline import *
from physics_toolkit import *
from physics_constants import *

import loadfile
import numpy

import chef_propagate

def read_em(index,units):
    x = loadfile.loadfile_transpose("code_parts_%d.dat" % index)
    shape = numpy.shape(x)
    for i in range(0,6):
        for j in range(0,shape[1]):
            x[i,j] *= units[i]
    return x

def compare(local_particles,index,units):
    leo = read_em(index,units)
    syn = local_particles
    passes = True
    sum = 0.0
    for i in range(0,6):
        for j in range(0,10):
            if leo[i,j] == 0.0:
                norm = 1.0
            else:
                norm = leo[i,j]
            norm_diff = (leo[i,j] - syn[i,j])/norm
            sum += norm_diff
            if (abs(norm_diff) > 1.0e-5):
                print i,j,norm_diff,leo[i,j],syn[i,j]
                passes = False
    print "average diff =", sum/60.0
    return passes

def read_parts(index,units):
    retval = numpy.zeros((7,10),'d')
    retval[0:6,:] = read_em(index,units)
    return retval

order=1
JetParticle.createStandardEnvironments(order)
factory = XSIFFactory("ertml_filecalls.xsif")
orig = factory.create_beamline("ERTML")
bmln = DriftsToSlots(orig)

energy = bmln.Energy()
energy = 5.0005110033999997654
positron = Positron ( energy )
rrv = RefRegVisitor ( positron )
bmln.accept( rrv )

units = numpy.array((0.1,0.2,0.3,0.4,0.5,-0.6),'d')
particles = read_parts(-1,units)
step = 0
for elem in bmln:
    chef_propagate.chef_propagate(particles,10,elem,
        energy,"positron",units,units)
    print step
    compare(particles,step,units)
    particles = read_parts(step,units)
    step += 1
