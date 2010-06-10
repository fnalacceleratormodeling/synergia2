#!/usr/bin/env python

from s2_fish import *
from macro_bunch import Macro_bunch
import numarray
import time
import pylab
import math
from math import pi
import sys

from test_solver import potential_uniform_sphere, compare_on_axis



def run():
    shape = (32,48,64)
    size = (2.0,5.0,3.0)
    offset = (0.1,0.5,0.4)

#     shape = (48,48,48)
#     size = (2.0,2.0,2.0)
#     offset = (0.0,0.0,0.0)

    sf = Real_scalar_field(shape,size,offset)
    mb = Macro_bunch()
    Q = 100000
    r0 = 0.2
    mb.init_sphere(Q,r0)
    total_charge = deposit_charge_cic(sf,mb.store)
    phi = solver_fft_open(sf)

    labels = ['x','y','z']
    for axis in range(0,3):
        r,phi_r,exact,maxerr,meanerr = compare_on_axis(axis,shape,size,offset,
                                                       phi,Q,r0)
        print labels[axis],": maxerr =",maxerr,"meanerr =",meanerr
        pylab.subplot(2,2,axis+1)
        pylab.plot(r,phi_r,'o')
        pylab.plot(r,exact,'r')
        pylab.xlabel(labels[axis])
    pylab.show()
run()
