#!/usr/bin/env python

from pytables_test_options import opts
import numpy
import tables

if __name__ == "__main__":
    print "hello world from tables_test.py"
    fileh = tables.openFile("array1.h5", mode = "w")
    root = fileh.root
    a = numpy.arange(24, dtype=numpy.float64).reshape(4,3,2)
    hdfarray = fileh.createArray(root, 'array_f', a, "3-D float array")
    fileh.close()
    
    fileh = tables.openFile("array1.h5", mode = "r")
    root = fileh.root
    
    newa = root.array_f.read()
    print "newa ="
    print newa
    fileh.close()
    print "pytables_test complete"