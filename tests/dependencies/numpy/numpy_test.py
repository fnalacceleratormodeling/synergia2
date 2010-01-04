#!/usr/bin/env python

from numpy_test_options import opts
import numpy

if __name__ == "__main__":
    print "hello world from numpy_test.py"
    z = numpy.zeros([10,2],'d')
    o = numpy.ones([5,4],'d')
    r = numpy.array(range(0,10),'d')
    print "zeros:"
    print z
    print "ones:"
    print o
    print "1-10:"
    print r
    print "numpy_test complete"


