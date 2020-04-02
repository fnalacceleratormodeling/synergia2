#!/usr/bin/env python

import sys
import os
import numpy as np

from nose.tools import *

testfiles = [
     "mpole_k1", "mpole_k1_tilt",
     "mpole_k2", "mpole_k2_tilt",
     "mpole_k3", "mpole_k3_tilt",
     "mpole_k4", "mpole_k4_tilt",
     "mpole_k5", "mpole_k5_tilt",
     "mpole_k6", "mpole_k6_tilt",
     "m_base_quad", "m_skew_quad", "m_tilt_quad"]

def test_file_outputs():
    for nm in testfiles:
        print("checking ", nm)
        m8p = np.load("m8"+nm+".npy")
        mxp = np.load(nm+".npy")
        assert(m8p.shape[0] == 32)
        assert(mxp.shape[0] == 32)
        numpart = mxp.shape[0]
        for p in range(numpart):
            for j in range(4):
                print("    particle:  [%d, %d] %.14g <-> %.14g"%(p,j, m8p[p,j], mxp[p,j]))
                assert_almost_equal(m8p[p, j], mxp[p, j], 12)
