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
    def test_01_get_cylindrical_coords(self):
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
        
    def xtest_03_cylindrical_field_domain_construct(self):
        radius = 2.0
        length = 1.5
        grid_shape = [32,32,32]
        periodic_z = True
        field_domain = Cylindrical_field_domain(radius,length,grid_shape,periodic_z)
        
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

    def test_07_aaarandom(self):
        print "jfa: duh"
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        num = 100
        radius = 5.0
        p = Numeric.zeros((7,num),'d')
        
        for i in range(0,num):
            for j in range(1,7):
                p[j,i] = 0.0
            p[0,i] = (radius*i)/num
        
        mb.units = Numeric.ones((6),'d')
        mb.local_num = num
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = p
        
        #mb.write_particles_text("cylinder")
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        grid_shape = [9,4,4]
        z_length = 2*pi
        field_domain = Cylindrical_field_domain(radius,z_length,grid_shape,True)
        rho = Numeric.zeros(grid_shape,'d')
        print "jfa: about to deposit"
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        sys.exit(999)
        print "jfa: deposited"
        nonzero = 0
        integral = 0
        cell_size = field_domain.get_cell_size()
        rhoz = Numeric.zeros((grid_shape[0],),'d')
        for r in range(0,grid_shape[0]):
            for theta in range(0,grid_shape[1]):
                for z in range(0,grid_shape[2]):
                    rhoz[r] += rho[r,theta,z]
                    if rho[r,theta,z] != 0:
                        nonzero += 1
                        integral += rho[r,theta,z]*0.5*((r+1)**2 - r**2)* \
                            (cell_size[0]**2*cell_size[1]*cell_size[2])
        print "found",nonzero,"entries. integral = %0.6f" % integral
        pylab.plot(rhoz)
        pylab.show()
            
    def xtest_07_aarandom(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        num = 100000
        p = Numeric.zeros((7,num),'d')
        covs = Numeric.zeros((6,6),'d')
        means = Numeric.zeros((6,),'d')

        for i in range(0,6):
            means[i] = 0.0;
            for j in range(0,6):
                if i==j:
                    covs[i,j] = 1.5

        t0 = time.time()
        populate.populate_uniform_cylinder(p,means,covs,0,0,False)
        t1 = time.time()

        mb.units = Numeric.ones((6),'d')
        mb.local_num = num
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = p
        
        #mb.write_particles_text("cylinder")
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        grid_shape = [9,4,4]
        radius = 5.0
        z_length = 2*pi
        field_domain = Cylindrical_field_domain(radius,z_length,grid_shape,True)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        nonzero = 0
        integral = 0
        cell_size = field_domain.get_cell_size()
        for r in range(0,grid_shape[0]):
            for theta in range(0,grid_shape[1]):
                for z in range(0,grid_shape[2]):
                    if rho[r,theta,z] != 0:
                        nonzero += 1
                        integral += rho[r,theta,z]*0.5*((r+1)**2 - r**2)* \
                            (cell_size[0]**2*cell_size[1]*cell_size[2])
        print "found",nonzero,"entries. integral = %0.6f" % integral
        for i in range(0,grid_shape[1]):
            for j in range(0,grid_shape[2]):
                pylab.plot(rho[:,i,j])
        pylab.show()
            
    def test_07_deposit(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        # jfa: this is a workaround for the lack of a reasonable general populate
        numcircs = 41
        numdisks = 21
        numtheta0 = 7
        num = 0
        for i in range(0,numcircs):
            num += numdisks*numtheta0*(i+1)
        print "num particles =",num
        p = Numeric.zeros((7,num),'d')

        mb.units = Numeric.ones((6),'d')
        mb.local_num = num
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = Numeric.zeros((7,mb.local_num),'d')
        populate.populate_uniform_cylinder_regular(mb.particles,2.0,2*pi,numcircs,numdisks,numtheta0)
        
        #mb.write_particles_text("cylinder")
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        grid_shape = [12,4,4]
        field_domain = Cylindrical_field_domain(5.0,2*pi,grid_shape,True)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        nonzero = 0
        integral = 0
        cell_size = field_domain.get_cell_size()
        for r in range(0,grid_shape[0]):
            for theta in range(0,grid_shape[1]):
                for z in range(0,grid_shape[2]):
                    if rho[r,theta,z] != 0:
                        nonzero += 1
                        integral += rho[r,theta,z]*0.5*((r+1)**2 - r**2)* \
                            (cell_size[0]**2*cell_size[1]*cell_size[2])
        print "found",nonzero,"entries. integral = %0.6f" % integral

        for i in range(0,grid_shape[1]):
            for j in range(0,grid_shape[2]):
                pylab.plot(rho[:,i,j])
        pylab.show()

    def xtest_07_deposit_integral(self):
        mb = Macro_bunch(physics_constants.PH_NORM_mp,1)
        # jfa: this is a workaround for the lack of a reasonable general populate
        num = 1
        p = Numeric.zeros((7,num),'d')

        mb.units = Numeric.ones((6),'d')
        mb.local_num = num
        mb.total_num = mb.local_num
        mb.ref_particle = Numeric.zeros((6,),'d')
        mb.ref_particle[5] = -1.1
        mb.is_fixed_z=0
        mb.total_current=1.0
        mb.charge=1

        mb.particles = Numeric.zeros((7,mb.local_num),'d')
        if len(sys.argv) > 1:
            x = float(sys.argv[1])
        else:
            x = -3.0
        if len(sys.argv) > 2:
            y = float(sys.argv[2])
        else:
            y = 1.0
        if len(sys.argv) > 3:
            z = float(sys.argv[3])
        else:
            z = 0.1
        print
        print "point =",x,y,z
        mb.particles[0,0] = x
        mb.particles[2,0] = y
        mb.particles[4,0] = z
        print "particles =",mb.particles
        #populate.populate_uniform_cylinder_regular(mb.particles,2.0,2*pi,numcircs,numdisks,numtheta0)
        
        #mb.write_particles_text("cylinder")
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        print "coords =",coords
        grid_shape = [8,6,6]
        radius = 5.0
        z_length = 2*pi
        field_domain = Cylindrical_field_domain(radius,z_length,grid_shape,True)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        #print "rho =",rho
        nonzero = 0
        integral = 0
        cell_size = field_domain.get_cell_size()
        print "cell_size =",cell_size
        for r in range(0,grid_shape[0]):
            for theta in range(0,grid_shape[1]):
                for z in range(0,grid_shape[2]):
                    if rho[r,theta,z] != 0:
                        nonzero += 1
                        integral += rho[r,theta,z]*0.5*((r+1)**2 - r**2)* \
                            (cell_size[0]**2*cell_size[1]*cell_size[2])
                        print "nonzero:",r,theta,z,":",rho[r,theta,z]
        print "found",nonzero,"entries. integral = %0.6f" % integral
        
    def xtest_08_tridiag(self):
        A = Numeric.array([[1,2,0,0],[3,4,5,0],[0,6,7,8],[0,0,9,10]],'D')
        # next three lines are the three tridiagonal vectors of A
        diag = Numeric.array([1,4,7,10],'D')
        above_diag = Numeric.array([2,5,8],'D')
        below_diag = Numeric.array([3,6,9],'D')
        
        b = Numeric.array([1,2,3,4],'D')
        x = Numeric.zeros([4],'D')
        solve_tridiag_nonsym(diag,above_diag,below_diag,b,x)
        newb = Numeric.matrixmultiply(A,x)
        for i in range(0,4):
            self.assertAlmostEqual(b[i].real,newb[i].real)
            self.assertAlmostEqual(b[i].imag,newb[i].imag)

    def xtest_09_solve(self):
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
        r0 = 0.4
        for i in range(0,6):
            covs[i,i] = r0**2/4
        #~ populate.populate_uniform_cylinder_quasi(mb.particles,means,covs,0)
        populate.populate_uniform_cylinder(mb.particles,means,covs,0,0,1)
        mb.write_particles("debug.h5")
        coords = Numeric.zeros((3,mb.local_num),'d')
        get_cylindrical_coords(mb.get_store(),coords)
        physical_size = [1.0,2*pi,2*pi]
        physical_offset = [physical_size[0]/2.0,physical_size[1]/2.0,0.0]
        grid_shape = [6,4,4]
        periodic = [False,True,True]
        field_domain = Field_domain(physical_size,physical_offset,grid_shape,periodic)
        rho = Numeric.zeros(grid_shape,'d')
        deposit_charge_cic_cylindrical(field_domain, rho ,mb.get_store(),coords)
        phi = Numeric.zeros(grid_shape,'d')
        solve_cylindrical_finite_periodic(field_domain,rho,phi)
        #~ print "phi =",phi
        #for mphi in range(0,grid_shape[1]):
            #pylab.plot(phi[:,mphi,0])
        #pylab.show()
    

if __name__ == '__main__':
    unsuccessful = 0
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_cylindrical)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
