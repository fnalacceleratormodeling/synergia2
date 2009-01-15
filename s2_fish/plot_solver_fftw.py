#!/usr/bin/env python

from s2_solver_fftw import *
from s2_containers import *
from s2_deposit import *
from s2_electric_field import *

import macro_bunch
import numarray
import time
import pylab
import math
from math import pi
import sys
import sublocal_paths
import physics_constants
import loadfile

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
#~ print "shape =",shape
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
print "about to solve on",MPI.COMM_WORLD.Get_rank()
fftwh = Fftw_helper(shape)
phi = solver_fftw_open(rho,fftwh,0)
print "completed solve on",MPI.COMM_WORLD.Get_rank()
sys.stdout.flush()
#~ phi.write_to_file("phi-%d-%d" % (MPI.COMM_WORLD.Get_size(),MPI.COMM_WORLD.Get_rank()))
#phi.get_points().print_("phi")
t_solve = time.time() - t0
print "about to plot on",MPI.COMM_WORLD.Get_rank()
sys.stdout.flush()
if MPI.COMM_WORLD.Get_rank()<MPI.COMM_WORLD.Get_size()/2 or MPI.COMM_WORLD.Get_size() == 1:
    #~ print "init =",t_init,", deposit =",t_deposit,", solve =",t_solve
    points = phi.get_points()
    shape = points.get_shape()
    print "on",MPI.COMM_WORLD.Get_rank(),"shape is",shape,points.get_dim0_lower(),points.get_dim0_upper()
    p = []
    x = []
    for i in range(points.get_dim0_lower(),points.get_dim0_upper()):
        x.append(float(i))
        index = (i,shape[1]/2,shape[2]/2)
        p.append(points.get(index))
    #~ if MPI.COMM_WORLD.Get_size() == 1:
        #~ loadfile.savevector("answer",p)
    #~ answer = loadfile.loadvector("answer")
    #~ pylab.plot(answer,'ro')
    #~ pylab.plot(x,p,'ro-')
    #~ analytic = []
    #~ for xi in x:
        #~ analytic.append(exact(Q,r0,
    num_points = n
    
    ax = numarray.arrayrange(num_points)*size[0]/(num_points -1) - 0.5*size[0]
    aexact = numarray.zeros([num_points],numarray.Float)
    aphi = numarray.zeros([num_points],numarray.Float)
    index = 0;
    for i in range(0,len(ax)):
        x = ax[i]
        if abs(x - 0.5*size[0]) < 1.0e-10:
            x = 0.5*size[0] - 1.0e-10
        rvec = (x,0,0)
        try:
            aphi[i] = phi.get_val(rvec)
        except:
            pass
        aexact[i] = exact(Q,r0,rvec)
        #~ arho[i] = rho.get_val(rvec)
    ###    print x,aphi[i]/aexact[i]
    pylab.plot(ax,aphi,'o',label='solver using %s grid' % str(shape))
    pylab.plot(ax,aexact,'r',label='analytic solution')

    pylab.xlabel('$x$')
    pylab.ylabel('$\phi(x,0,0)$')
    pylab.legend(loc=0)
    sys.stdout.flush()
    if MPI.COMM_WORLD.Get_rank() < MPI.COMM_WORLD.Get_size()/2 or MPI.COMM_WORLD.Get_size() == 1:
        print "doing show on rank",MPI.COMM_WORLD.Get_rank()
        #~ pylab.savefig("foo%d.eps" % MPI.COMM_WORLD.Get_rank())
        pylab.show()
MPI.WORLD.Barrier()
print "ended successfully on",MPI.COMM_WORLD.Get_rank()

