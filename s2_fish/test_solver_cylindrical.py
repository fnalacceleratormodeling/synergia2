#!/usr/bin/env python


from synergia import physics_constants
from macro_bunch import Macro_bunch
from s2_solver_cylindrical import *
import populate

import Numeric
import MLab

import unittest
import sys
from math import cos,sin,pi

#~ begin debug
import pylab
import time
#~ end debug

class Test_solver_cylindrical(unittest.TestCase):                       
    def xtest_01_get_cylindrical_coords(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        Q = 10
        r0 = 0.2
        arbitrary_z = -1.9
        mb.init_sphere(Q,r0)
        for i in range(0,Q):
            theta = 6.0/Q*i
            mb.get_local_particles()[0,i] = r0*cos(theta)
            mb.get_local_particles()[2,i] = r0*sin(theta)
            mb.get_local_particles()[4,i] = arbitrary_z
        local_num = mb.get_num_particles_local()
        coords = Numeric.zeros((3,local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        for i in range(0,Q):
            theta = 6.0/Q*i
            self.assertAlmostEqual(coords[0,i],r0)
            self.assertAlmostEqual(coords[1,i],theta)
            self.assertAlmostEqual(coords[2,i],arbitrary_z)
    
    def xtest_02_field_domain_construct1(self):
        field_domain = Field_domain()
        physical_size = [2,4,6]
        physical_offset = [0.1,0.2,0.3]
        grid_shape = [32,32,32]
        periodic = [False,True,True]
        field_domain.set_params(physical_size,physical_offset,grid_shape,periodic)
    
    def xtest_03_field_domain_construct2(self):
        physical_size = [2,4,6]
        physical_offset = [0.1,0.2,0.3]
        grid_shape = [32,32,32]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        
    def xtest_04_field_domain_get_grid_shape(self):
        physical_size = [2,4,6]
        physical_offset = [0.1,0.2,0.3]
        grid_shape = [32,32,32]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        returned_shape = field_domain.get_grid_shape()
        for i in range(0,3):
            self.assertEqual(returned_shape[i],grid_shape[i])

    def xtest_05_field_domain_get_cell_size(self):
        physical_size = [2,4,6]
        physical_offset = [0.1,0.2,0.3]
        grid_shape = [32,32,32]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        cell_size= field_domain.get_cell_size()
        for i in range(0,3):
            expected = physical_size[i]/(grid_shape[i] - 1.0)
            self.assertAlmostEqual(cell_size[i],expected)

    def xtest_06_field_domain_get_periodic(self):
        physical_size = [2,4,6]
        physical_offset = [0.1,0.2,0.3]
        grid_shape = [32,32,32]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        returned_periodic = field_domain.get_periodic()
        for i in range(0,3):
            self.assertEqual(returned_periodic[i],periodic[i])

    def xtest_07_deposit(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        # jfa: this is a workaround for the lack of a reasonable general populate
        mb.units = Numeric.ones((6),'d')
        mb.local_num = 20000
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = Numeric.zeros((7,mb.local_num),'d')
        means = Numeric.zeros((6),'d')
        covs = Numeric.zeros((6,6),'d')
        for i in range(0,6):
            covs[i,i] = 1.0
        #~ populate.populate_uniform_cylinder_quasi(mb.particles,means,covs,0)
        populate.populate_uniform_cylinder(mb.particles,means,covs,0,0,1)
        
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        physical_size = [5.0,2*pi,2*pi]
        physical_offset = [physical_size[0]/2.0,physical_size[1]/2.0,0.0]
        grid_shape = [20,8,4]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        # jfa: need to test resulting rho...
        
    def xtest_08_solve(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        # jfa: this is a workaround for the lack of a reasonable general populate
        mb.units = Numeric.ones((6),'d')
        mb.local_num = 20000
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = Numeric.zeros((7,mb.local_num),'d')
        means = Numeric.zeros((6),'d')
        covs = Numeric.zeros((6,6),'d')
        for i in range(0,6):
            covs[i,i] = 1.0
        #~ populate.populate_uniform_cylinder_quasi(mb.particles,means,covs,0)
        populate.populate_uniform_cylinder(mb.particles,means,covs,0,0,1)
        
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        physical_size = [5.0,2*pi,2*pi]
        physical_offset = [physical_size[0]/2.0,physical_size[1]/2.0,0.0]
        grid_shape = [20,8,4]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        phi = Numeric.zeros(grid_shape,'d')
        solve_cylindrical_finite_periodic(field_domain,rho,phi)

if __name__ == '__main__':
    unsuccessful = 0
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_cylindrical)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
