#!/usr/bin/env python

import populate
import Numeric
import MLab
from ArrayPrinter import array2string

import unittest
import sys
import time

import pylab

class Test_populate(unittest.TestCase):
    def test_01_populate_6d_gaussian(self):
        N = 10000
        p = Numeric.zeros((7,N),'d')
        covs = Numeric.zeros((6,6),'d')
        means = Numeric.zeros((6,),'d')

        for i in range(0,6):
            means[i] = 0.01*i;
            for j in range(0,6):
                if i==j:
                    covs[i,j] = i+1.0
                if (i+1 == j) and (i%2 == 0):
                    covs[i,j] = 0.1*(j+1.0)
                    covs[j,i] = 0.1*(j+1.0)

        #~ print "requested covs"
        #~ print covs
        #~ print "requested means"
        #~ print means
        t0 = time.time()
        populate.populate_6d_gaussian(p,means,covs,0,0,False)
        t1 = time.time()
        print "\npopulate took",t1-t0

        # covs and means were overwritten by populate_6d_gaussian
        for i in range(0,6):
            means[i] = 0.01*i;
            for j in range(0,6):
                if i==j:
                    covs[i,j] = i+1.0
                if (i+1 == j) and (i%2 == 0):
                    covs[i,j] = 0.1*(j+1.0)
                    covs[j,i] = 0.1*(j+1.0)

        found_covs = MLab.cov(Numeric.transpose(p))
        found_means = MLab.mean(Numeric.transpose(p))
        
        #~ print "actual covariances"
        #~ print array2string(found_covs[0:6,0:6],suppress_small=1,precision=4)
        #~ print "actual means"
        #~ print found_means

        for i in range(0,6):
            self.assertAlmostEqual(means[i],found_means[i])
            for j in range(0,6):
                self.assertAlmostEqual(covs[i,j],found_covs[i,j]*(N-1.0)/N) # N vs. N-1!!!!!!!
        
        pylab.figure(0)
        pylab.plot(p[0,:],p[2,:],'.')

    def test_02_populate_6d_gaussian_quasi(self):
        N = 10000
        p = Numeric.zeros((7,N),'d')
        covs = Numeric.zeros((6,6),'d')
        means = Numeric.zeros((6,),'d')

        for i in range(0,6):
            means[i] = 0.01*i;
            for j in range(0,6):
                if i==j:
                    covs[i,j] = i+1.0
                if (i+1 == j) and (i%2 == 0):
                    covs[i,j] = 0.1*(j+1.0)
                    covs[j,i] = 0.1*(j+1.0)

        #~ print "requested covs"
        #~ print covs
        #~ print "requested means"
        #~ print means
        t0 = time.time()
        populate.populate_6d_gaussian_quasi(p,means,covs,0)
        t1 = time.time()
        print "\nquasi populate took",t1-t0

        # covs and means were overwritten by populate_6d_gaussian
        for i in range(0,6):
            means[i] = 0.01*i;
            for j in range(0,6):
                if i==j:
                    covs[i,j] = i+1.0
                if (i+1 == j) and (i%2 == 0):
                    covs[i,j] = 0.1*(j+1.0)
                    covs[j,i] = 0.1*(j+1.0)

        found_covs = MLab.cov(Numeric.transpose(p))
        found_means = MLab.mean(Numeric.transpose(p))
        
        #~ print "actual covariances"
        #~ print array2string(found_covs[0:6,0:6],suppress_small=1,precision=4)
        #~ print "actual means"
        #~ print found_means

        for i in range(0,6):
            self.assertAlmostEqual(means[i],found_means[i])
            for j in range(0,6):
                self.assertAlmostEqual(covs[i,j],found_covs[i,j]*(N-1.0)/N) # N vs. N-1!!!!!!!

        pylab.figure(1)
        pylab.plot(p[0,:],p[2,:],'.')
        pylab.show()

if __name__ == '__main__':
    unsuccessful = 0
    populate_suite = unittest.TestLoader().loadTestsFromTestCase(Test_populate)
    retval = unittest.TextTestRunner(verbosity=2).run(populate_suite)
    if not retval.wasSuccessful():
        unsuccessful = 1

    sys.exit(unsuccessful)
