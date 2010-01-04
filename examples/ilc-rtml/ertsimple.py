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

def read_em(particles,index):
    x = loadfile.loadfile_transpose("code_parts_%d.dat" % index)
    shape = numpy.shape(x)
    for p in range(0,len(particles)):
        particles[p].set_x(x[0,p])
        particles[p].set_npx(x[1,p])
        particles[p].set_y(x[2,p])
        particles[p].set_npy(x[3,p])
        particles[p].set_cdt(x[4,p])
        particles[p].set_ndp(x[5,p])

def parts_to_array(particles):
    retval = numpy.zeros((6,len(particles)),'d')
    for p in range(0,len(particles)):
        s = particles[p]
        retval[0,p] = s.get_x()
        retval[1,p] = s.get_npx()
        retval[2,p] = s.get_y()
        retval[3,p] = s.get_npy()
        retval[4,p] = s.get_cdt()
        retval[5,p] = s.get_ndp()
    return retval
    
def compare(particles,index):
    leo = loadfile.loadfile_transpose("code_parts_%d.dat" % index)
    syn = parts_to_array(particles)
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
    
factory = XSIFFactory("ertml_filecalls.xsif")
orig = factory.create_beamline("ERTML")
bmln = DriftsToSlots(orig)

energy = bmln.Energy()
positron = Positron ( energy )
rrv = RefRegVisitor ( positron )
bmln.accept( rrv )

numpart = 10
particles = []
for i in range(0,numpart):
    particles.append(Positron(energy))

read_em(particles,-1)
#~ for p in particles:
    #~ print p.State()

step = 0
for elem in bmln:
    part = 0
    for particle in particles:
        print "before part",part,":",particle.State()
        elem.propagateParticle(particle)
        print "after part",part,":",particle.State()
        part += 1
    print step
    compare(particles,step)
    read_em(particles,step)
    step += 1
