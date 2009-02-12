#!/usr/bin/env python

from s2_solver_fftw3 import *
from s2_containers import *
from s2_deposit import *
from s2_electric_field import *

import macro_bunch
import numpy
import time
import pylab
import math
from math import pi
import sys


def exact(Q,r0,rvec):
    r = math.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    if r <= r0:
        retval = Q/(8*pi*r0**3)*(3*r0**2-r**2)
    else:
        retval = Q/(4*pi*r)
    return retval

n = int(sys.argv[1])
shape = (n,n,n)
print "shape =",shape
size = (2.0,2.0,2.0)
#size = (100.0,100.0,100.0)
rho = Real_scalar_field(shape,size,(0.0,0.0,0.0))
mb = macro_bunch.Macro_bunch()
t0 = time.time()
Q = 100000
r0 = 0.5
mb.init_sphere(Q,r0)
t_init = time.time() - t0
t0 = time.time()
total_charge = deposit_charge_cic(rho,mb.get_store(),0)
t_deposit = time.time() - t0
t0 = time.time()
phi = solver_fftw3_open(rho,0)
#phi.get_points().print_("phi")
t_solve = time.time() - t0

print "init =",t_init,", deposit =",t_deposit,", solve =",t_solve
num_points = n
ax = numpy.arrayrange(num_points)*size[0]/(num_points -1) - 0.5*size[0]
aphi = numpy.zeros([num_points],numpy.Float)
arho = numpy.zeros([num_points],numpy.Float)
aexact = numpy.zeros([num_points],numpy.Float)
index = 0;
for i in range(0,len(ax)):
    x = ax[i]
    if abs(x - 0.5*size[0]) < 1.0e-10:
        x = 0.5*size[0] - 1.0e-10
    rvec = (x,0,0)
    aphi[i] = phi.get_val(rvec)
    aexact[i] = exact(Q,r0,rvec)
    arho[i] = rho.get_val(rvec)
###    print x,aphi[i]/aexact[i]
pylab.plot(ax,aphi,'o')
pylab.plot(ax,aexact,'r')

# ax2 = pylab.twinx()
# pylab.plot(ax,aexact/aphi,'g')
# ax2.yaxis.tick_right()

# pylab.figure()
# pylab.plot(ax,arho)
#print aphi[0]/aexact[0], aphi[0]/aexact[0]*4*pi, aphi[0]/aexact[0]/(pi)*n*n/(4*32768.0)
print aphi[0]/aexact[0], aphi[0]/aexact[0]*4*pi, aphi[0]/aexact[0]
sys.stdout.flush()
pylab.show()

