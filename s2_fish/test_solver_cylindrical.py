#!/usr/bin/env python

from s2_solver_cylindrical import *

from synergia import physics_constants
from macro_bunch import Macro_bunch

import Numeric

import unittest
import sys
from math import cos,sin

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
                   
if __name__ == '__main__':
    unsuccessful = 0
    solver_suite = unittest.TestLoader().loadTestsFromTestCase(Test_solver_cylindrical)
    retval = unittest.TextTestRunner(verbosity=2).run(solver_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
