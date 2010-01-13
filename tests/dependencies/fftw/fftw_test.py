#!/usr/bin/env python

from fftw_test_options import opts
import trivial_fftw

if __name__ == "__main__":
    print "about to call C++ function:"
    trivial_fftw.doit()
    print "fftw_test.py finished successfully"
