#!/usr/bin/env python

import s2_fish
import macro_bunch
import numarray
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

periodic = 1
n = int(sys.argv[1])
shape = (n,n,n)
print "shape =",shape
size = (2.0,2.0,2.0)
#size = (100.0,100.0,100.0)
rho = s2_fish.Real_scalar_field(shape,size,(0.0,0.0,0.0))
mb = macro_bunch.Macro_bunch()
t0 = time.time()
Q = 100000
r0 = 0.5
mb.init_cylinder(Q,r0,size[2])
t_init = time.time() - t0
t0 = time.time()
total_charge = s2_fish.deposit_charge_cic(rho,mb.get_store(),periodic)
print "total_charge =",total_charge
t_deposit = time.time() - t0
t0 = time.time()
phi = s2_fish.solver_fft_open(rho,periodic)
#phi.get_points().print_("phi")
t_solve = time.time() - t0

print "init =",t_init,", deposit =",t_deposit,", solve =",t_solve
num_points = n
ax = numarray.arrayrange(num_points)*size[0]/(num_points -1) - 0.5*size[0]
az = numarray.arrayrange(num_points)*size[2]/(num_points -1) - 0.5*size[2]
aphi = numarray.zeros([num_points],numarray.Float)
arho = numarray.zeros([num_points],numarray.Float)
aexact = numarray.zeros([num_points],numarray.Float)
index = 0;
for k in range(0,len(ax)):
    z = az[k]
    if abs(z - 0.5*size[2]) < 1.0e-10:
        z = 0.5*size[2] - 1.0e-10
    for i in range(0,len(ax)):
        x = ax[i]
        if abs(x - 0.5*size[0]) < 1.0e-10:
            x = 0.5*size[0] - 1.0e-10
        rvec = (x,0,z)
        aphi[i] = phi.get_val(rvec)
        aexact[i] = exact(Q,r0,rvec)
        arho[i] = rho.get_val(rvec)
###        print x,aphi[i]/aexact[i]
#    print "k=",k,aphi
    pylab.plot(ax,aphi)
###pylab.plot(ax,aexact,'r')

# ax2 = pylab.twinx()
# pylab.plot(ax,aexact/aphi,'g')
# ax2.yaxis.tick_right()

# pylab.figure()
# pylab.plot(ax,arho)
#print aphi[0]/aexact[0], aphi[0]/aexact[0]*4*pi, aphi[0]/aexact[0]/(pi)*n*n/(4*32768.0)
sys.stdout.flush()
pylab.show()

