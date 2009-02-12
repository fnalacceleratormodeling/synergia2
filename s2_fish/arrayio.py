#!/usr/bin/env python

import sys
import numpy

def readfort(filename,shape):
    expected_length = 1
    for si in shape:
        expected_length*=si
    data = []
    f = open(filename,"r")
    for line in f.readlines():
        data+=(map(float,line.split()))
    f.close()
    if not len(data) == expected_length:
        raise RuntimeError,"found length %d,\nbut expected an array of shape %s (length %d)" % (len(data),str(shape),expected_length)
    retval = numpy.zeros(shape,numpy.Float)
    index = 0
    for k in range(0,shape[2]):
        for j in range(0,shape[1]):
            for i in range(0,shape[0]):
                retval[i,j,k] = data[index]
                index += 1
    return retval

def readtensor(filename):
    data = []
    f = open(filename,"r")
    shape = map(int,f.readline().split())
    for line in f.readlines():
        data+=(map(float,line.split()))
    f.close()
    retval = numpy.zeros(shape,numpy.Float)
    index = 0
    for i in range(0,shape[0]):
        for j in range(0,shape[1]):
            for k in range(0,shape[2]):
                retval[i,j,k] = data[index]
                index += 1
    return retval

def readscalarfield(filename):
    data = []
    f = open(filename,"r")
    f.readline() # skip over scalar field information
    f.readline() # skip over scalar field information
    shape = map(int,f.readline().split())
    for line in f.readlines():
        data+=(map(float,line.split()))
    f.close()
    retval = numpy.zeros(shape,numpy.Float)
    index = 0
    for i in range(0,shape[0]):
        for j in range(0,shape[1]):
            for k in range(0,shape[2]):
                retval[i,j,k] = data[index]
                index += 1
    return retval

def readpetscvec(filename):
    data = []
    f = open(filename,"r")
    shape = map(int,f.readline().split())
    for line in f.readlines():
        cdata=line.split()
        for cdatum in cdata:
            pair = cdatum[1:len(cdatum)-1].split(',')
            data.append(complex(float(pair[0]),float(pair[1])))
    f.close()
    retval = numpy.zeros(shape,numpy.Complex)
    index = 0
    for i in range(0,shape[0]):
        for j in range(0,shape[1]):
            for k in range(0,shape[2]):
                retval[i,j,k] = data[index]
                index += 1
    return retval

if __name__ == "__main__":
#    foo = readfort("../fort.40",[8,8,8])
#    foo = readscalarfield("rho")
    foo = readpetscvec("G_hat2_petsc")
    print foo
    
