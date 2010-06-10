#!/usr/bin/env python

from s2_solver_fftw import *
from s2_containers import *
from s2_deposit import *
from s2_electric_field import *

import sublocal_paths
import physics_constants

import macro_bunch
import numarray
import time
import pylab
import math
from math import pi
import sys

from mpi4py import MPI

def exact(Q,r0,rvec):
    r = math.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    if r <= r0:
        retval = Q/(8*pi*r0**3)*(3*r0**2-r**2)
    else:
        retval = Q/(4*pi*r)
    return retval

n = int(sys.argv[1])
shape = (n,n,n)
size = (2.0,2.0,2.0)
#size = (100.0,100.0,100.0)
rho = Real_scalar_field(shape,size,(0.0,0.0,0.0))
mb = macro_bunch.Macro_bunch(physics_constants.PH_NORM_mp,1)
t0 = time.time()
Q = 100000
r0 = 0.5
mb.init_sphere(Q,r0)
t_init = time.time() - t0
t0 = time.time()
total_charge = deposit_charge_cic(rho,mb.get_store(),0)
#~ if MPI.COMM_WORLD.Get_rank() == 0:
    #~ for i in range(0,2):
        #~ for j in range(0,2):
            #~ for k in range(0,2):
                #~ rho.get_points().set((i,j,k),1.0)
                #~ print i,j,k,rho.get_points().get((i,j,k))
t_deposit = time.time() - t0
t0 = time.time()
fftwh = Fftw_helper(shape)
phi = solver_fftw_open(rho,fftwh,0)
phi.write_to_file("phi-%d-%d" % (MPI.COMM_WORLD.Get_size(),MPI.COMM_WORLD.Get_rank()))
#phi.get_points().print_("phi")
t_solve = time.time() - t0
E = calculate_E_n(phi,0)

#~ E = phi
if MPI.COMM_WORLD.Get_rank()<4:
    print "init =",t_init,", deposit =",t_deposit,", solve =",t_solve
    points = E.get_points()
    shape = points.get_shape()
    p = []
    for i in range(points.get_dim0_lower(),points.get_dim0_upper()):
        index = (i,shape[1]/2,shape[2]/2)
        p.append(points.get(index))
    pylab.plot(p,'o')
    pylab.title('proc %d' % MPI.COMM_WORLD.Get_rank())
    #~ num_points = n
    #~ ax = numarray.arrayrange(num_points)*size[0]/(num_points -1) - 0.5*size[0]
    #~ aphi = numarray.zeros([num_points],numarray.Float)
    #~ arho = numarray.zeros([num_points],numarray.Float)
    #~ aexact = numarray.zeros([num_points],numarray.Float)
    #~ index = 0;
    #~ for i in range(0,len(ax)):
        #~ x = ax[i]
        #~ if abs(x - 0.5*size[0]) < 1.0e-10:
            #~ x = 0.5*size[0] - 1.0e-10
        #~ rvec = (x,0,0)
        #~ try:
            #~ aphi[i] = phi.get_val(rvec)
        #~ except:
            #~ pass
        #~ aexact[i] = exact(Q,r0,rvec)
        #~ arho[i] = rho.get_val(rvec)
    #~ ###    print x,aphi[i]/aexact[i]
    #~ pylab.plot(ax,aphi,'o')
    #~ pylab.plot(ax,aexact,'r')


    sys.stdout.flush()
    if MPI.COMM_WORLD.Get_rank() < 4:
        print "doing show on rank",MPI.COMM_WORLD.Get_rank()
        pylab.show()

