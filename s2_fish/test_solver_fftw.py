#!/usr/bin/env python

from s2_solver_fftw import *
from s2_containers import *
from s2_deposit import *
from s2_electric_field import *

from synergia import physics_constants
from macro_bunch import Macro_bunch
import numarray
import time
import sys
import math
from math import pi,sin,cos

import unittest

def potential_uniform_sphere(Q,r0,rvec):
    """ The potential due to a uniform sphere of radius r0 and charge
    Q, centered at the origin."""
    r = math.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    if r <= r0:
        retval = Q/(8*pi*r0**3)*(3*r0**2-r**2)
    else:
        retval = Q/(4*pi*r)
    return retval

def potential_uniform_sphere_periodic(Q,r0,rvec,period):
    retval = 0.0
    num_images = 20
    image_vec = [0.0,0.0,0.0]
    image_vec[0] = rvec[0]
    image_vec[1] = rvec[1]
    for image in range(-num_images,num_images+1):
        image_vec[2] = rvec[2] - image*period
        retval += potential_uniform_sphere(Q,r0,image_vec)
    return retval

def E_field_uniform_sphere(Q,r0,rvec,axis):
    """ The electric field in the axis direction due to a uniform
    sphere of radius r0 and charge Q, centered at the origin."""
    r = math.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    rhat = [rvec[0]/r,rvec[1]/r,rvec[2]/r]
    if r <= r0:
        dphi_dr = Q/(8*pi*r0**3)*(-2*r)
    else:
        dphi_dr = -Q/(4*pi*r**2)
    return dphi_dr*rhat[axis]

def E_field_uniform_sphere_periodic(Q,r0,rvec,period,axis):
    retval = 0.0
    num_images = 10
    image_vec = [0.0,0.0,0.0]
    image_vec[0] = rvec[0]
    image_vec[1] = rvec[1]
    for image in range(-num_images,num_images+1):
        image_vec[2] = rvec[2] - image*period
        retval += E_field_uniform_sphere(Q,r0,image_vec,axis)
    return retval

def compare_on_axis(axis,shape,size,offset,phi,Q,r0,periodic=False):
    r = numarray.arrayrange(shape[axis])*size[axis]/(shape[axis]-1) \
        - 0.5*size[axis] + offset[axis]
    phi_r = numarray.zeros(shape[axis],numarray.Float)
    exact = numarray.zeros(shape[axis],numarray.Float)
    max_err = 0.0
    sum_err = 0.0
    rfile = open("r.dat","w")
    phifile = open("phi.dat","w")
    exactfile = open("exact.dat","w")
    for i in range(0,shape[axis]):
        rvec = [0.0,0.0,0.0]
        ri = r[i]
        if abs(ri - 0.5*size[axis] - offset[axis]) < 1.0e-10:
            ri = 0.5*size[axis] + offset[axis]- 1.0e-10
        rvec[axis] = ri
        phi_r[i] = phi.get_val(rvec)
        if periodic:
            exact[i] = potential_uniform_sphere_periodic(Q,r0,rvec,size[2])
        else:
            exact[i] = potential_uniform_sphere(Q,r0,rvec)
        abs_err = abs(1.0 - phi_r[i]/exact[i])
        if abs_err > max_err:
            max_err = abs_err
        sum_err += abs_err
        rfile.write("%g\n" % r[i])
        phifile.write("%g\n" % phi_r[i])
        exactfile.write("%g\n" % exact[i])
    rfile.close()
    phifile.close()
    exactfile.close()
    mean_err = sum_err/shape[axis]
    return r,phi_r,exact,max_err,mean_err

def compare_E_on_axis(axis,shape,size,offset,E,Q,r0,E_axis,periodic=False):
    r = numarray.arrayrange(shape[axis])*size[axis]/(shape[axis]-1) \
        - 0.5*size[axis] + offset[axis]
    Er = numarray.zeros(shape[axis],numarray.Float)
    exact = numarray.zeros(shape[axis],numarray.Float)
    max_err = 0.0
    sum_err = 0.0
    rfile = open("er.dat","w")
    efile = open("e.dat","w")
    exactfile = open("eexact.dat","w")
    for i in range(0,shape[axis]):
        rvec = [0.0,0.0,0.0]
        ri = r[i]
        if abs(ri - 0.5*size[axis] - offset[axis]) < 1.0e-10:
            ri = 0.5*size[axis] + offset[axis]- 1.0e-10
        rvec[axis] = ri
        Er[i] = E.get_val(rvec)
        if periodic:
            exact[i] = E_field_uniform_sphere_periodic(Q,r0,rvec,size[2],E_axis)
        else:
            exact[i] = E_field_uniform_sphere(Q,r0,rvec,E_axis)
        if (exact[i] == 0.0):
            abs_err = abs(Er[i] - exact[i])
        else:
            abs_err = abs(1.0 - Er[i]/exact[i])
        if abs_err > max_err:
            max_err = abs_err
        sum_err += abs_err
        rfile.write("%g\n" % r[i])
        efile.write("%g\n" % Er[i])
        exactfile.write("%g\n" % exact[i])
    rfile.close()
    efile.close()
    exactfile.close()
    mean_err = sum_err/shape[axis]
    return r,Er,exact,max_err,mean_err

class Test_solver_fftw_open(unittest.TestCase):    
    def test_01_rough_grid(self):
        shape = (16,16,16)
        size = (2.0,2.0,2.0)
        offset = (0.0,0.0,0.0)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        fftwh = Fftw_helper(shape,False)
        phi = solver_fftw_open(sf,fftwh,0)
        for axis in range(0,3):
            r,phi_r,exact,max_err,mean_err = compare_on_axis(axis,shape,size,
                                                             offset,phi,Q,r0)
            self.failIf(max_err>0.07)
            self.failIf(mean_err>0.02)

    def test_02_fine_grid(self):
        shape = (48,48,48)
        size = (2.0,2.0,2.0)
        offset = (0.0,0.0,0.0)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        fftwh = Fftw_helper(shape,False)
        phi = solver_fftw_open(sf,fftwh,0)
        for axis in range(0,3):
            r,phi_r,exact,max_err,mean_err = compare_on_axis(axis,shape,size,
                                                             offset,phi,Q,r0)
            self.failIf(max_err>0.015)
            self.failIf(mean_err>0.003)

    def test_03_asymmetric_grid(self):
        shape = (16,24,32)
        size = (2.0,5.0,3.0)
        offset = (0.1,0.5,0.4)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        fftwh = Fftw_helper(shape,False)
        phi = solver_fftw_open(sf,fftwh,0)
        max_tolerance = [0.07,0.05,0.09]
        mean_tolerance = [0.03,0.006,0.03]
        for axis in range(0,3):
            r,phi_r,exact,max_err,mean_err = compare_on_axis(axis,shape,size,
                                                             offset,phi,Q,r0)
            self.failIf(max_err>max_tolerance[axis])
            self.failIf(mean_err>mean_tolerance[axis])

    def test_04_E_asymmetric_grid(self):
        shape = (24,32,48)
        size = (2.0,5.0,3.0)
        offset = (0.1,0.5,0.4)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        fftwh = Fftw_helper(shape,False)
        phi = solver_fftw_open(sf,fftwh,0)
        max_tolerance = [[0.3,   2.5e3, 2.5e3],
                         [4.5e3,   0.3, 4.5e3],
                         [2.5e3, 2.5e3,   0.3]]

        mean_tolerance = [[0.06,  500.0,  500.0],
                          [2.0e3,  0.08,  2.0e3],
                          [5000.0, 5000.0,   0.06]]
        for E_axis in range(0,3):
            E = calculate_E_n(phi,E_axis)
            for axis in range(0,3):
                r,Er,exact,max_err,mean_err = \
                                           compare_E_on_axis(axis,
                                                             shape,size,
                                                             offset,E,Q,
                                                             r0,E_axis)
                self.failIf(max_err>max_tolerance[E_axis][axis])
                self.failIf(mean_err>mean_tolerance[E_axis][axis])

    def test_05_kick(self):
        shape = (48,48,48)
        size = (2.0,2.0,2.0)
        offset = (0.0,0.0,0.0)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),0)
        fftwh = Fftw_helper(shape,False)
        phi = solver_fftw_open(sf,fftwh,0)
        for E_axis in range(0,3):
            E = calculate_E_n(phi,E_axis)
            apply_E_n_kick(E,E_axis,1.0,mb.get_store())

class Test_solver_fftw_open_periodic(unittest.TestCase):    
    #~ def test_01_rough_grid(self):
        #~ shape = (16,16,16)
        #~ size = (2.0,2.0,2.0)
        #~ offset = (0.0,0.0,0.0)
        #~ sf = Real_scalar_field(shape,size,offset)
        #~ mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        #~ Q = 100000
        #~ r0 = 0.2
        #~ mb.init_sphere(Q,r0)
        #~ total_charge = deposit_charge_cic(sf,mb.get_store(),1)
        #~ fftwh = Fftw_helper(shape,False)
        #~ phi = solver_fftw_open(sf,fftwh,1)
        #~ for axis in range(0,3):
            #~ r,phi_r,exact,max_err,mean_err = compare_on_axis(axis,shape,size,
                                                             #~ offset,phi,Q,r0,periodic=True)
            #~ self.failIf(max_err>0.07)
            #~ self.failIf(mean_err>0.02)

    def test_02_rough_grid_E(self):
        shape = (16,16,16)
        shape = (64,64,64)
        size = (2.0,2.0,2.0)
        offset = (0.0,0.0,0.0)
        sf = Real_scalar_field(shape,size,offset)
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 100000
        r0 = 0.2
        mb.init_sphere(Q,r0)
        total_charge = deposit_charge_cic(sf,mb.get_store(),1)
        fftwh = Fftw_helper(shape,True)
        phi = solver_fftw_open(sf,fftwh,1)
        for axis in range(0,3):
            E = calculate_E_n(phi,axis)
            r,E_r,exact,max_err,mean_err = compare_E_on_axis(axis,shape,size,
                                                             offset,E,Q,r0,periodic=True,E_axis=axis)
            print "axis:",axis,"max_err:",max_err,"mean_err:",mean_err
            #~ self.failIf(max_err>0.07)
            #~ self.failIf(mean_err>0.02)

if __name__ == '__main__':
    unsuccessful = 0
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_fftw_open)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1
    solver_suite_periodic = unittest.TestLoader().loadTestsFromTestCase(Test_solver_fftw_open_periodic)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite_periodic)
    if not retval.wasSuccessful():
        unsuccessful = 1

#     fd_solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_fd_multigrid)
#     retval = unittest.TextTestRunner(verbosity=2).run(fd_solver_suite)
#     if not retval.wasSuccessful():
#         unsuccessful = 1

    sys.exit(unsuccessful)
