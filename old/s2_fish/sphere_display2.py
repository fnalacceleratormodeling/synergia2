#!/usr/bin/env python

import pylab
import matplotlib
import math
import numarray
import sys

class Py_scalar_field:
    def __init__(self):
        self.physical_size = None
        self.physical_offset = None
        self.points = None

    def h(self):
        n = numarray.array(self.points.shape,numarray.Float)
        return self.physical_size/(n - 1.0)

    def axis(self,index):
        len = self.points.shape[index]
        retval = numarray.zeros([len],numarray.Float)
        for i in range(0,len):
            retval[i] = i*self.h()[index] + self.physical_offset[index] -\
                        self.physical_size[index]/2.0
        return retval
    
def readscalarfield(filename):
    sf = Py_scalar_field()
    data = []
    f = open(filename,"r")
    sf.physical_size = numarray.array(map(float,f.readline().split()))
    sf.physical_offset = numarray.array(map(float,f.readline().split()))
    shape = map(int,f.readline().split())
    for line in f.readlines():
        data+=(map(float,line.split()))
    f.close()
    sf.points = numarray.zeros(shape,numarray.Float)
    index = 0
    for i in range(0,shape[0]):
        for j in range(0,shape[1]):
            for k in range(0,shape[2]):
                sf.points[i,j,k] = data[index]
                index += 1
    return sf

def display161616(sf):
    index = 1
    max = sf.max()
    for i in range(0,4):
        for j in range(0,4):
            pylab.subplot(4,4,index)
            pylab.pcolor(sf[:,:,index-1],
                         norm=matplotlib.colors.normalize(vmax=max,vmin=0),
                         vmax=max,vmin=0.0)
            index +=1

def display323232(sf):
    index = 1
    max = sf.points.max()
    for i in range(0,32):
        pylab.subplot(6,6,index)
        pylab.pcolor(sf.points[:,:,index-1],
                     norm=matplotlib.colors.normalize(vmax=max,vmin=0),
                     vmax=max,vmin=0.0)
        index +=1

def display_axes(sf):
    shape = sf.points.shape;
    index = [shape[0]/2,shape[1]/2,shape[2]/2]
    pylab.subplot(2,2,1)
    pylab.plot(sf.axis(0),sf.points[:,index[1],index[2]],'.')
    pylab.subplot(2,2,2)
    pylab.plot(sf.axis(1),sf.points[index[0],:,index[2]],'.')
    pylab.subplot(2,2,3)
    pylab.plot(sf.axis(2),sf.points[index[0],index[1],:],'.')

def display_theory(sf):
    Q = 100000
    R = 0.1
    for i in range(0,3):
        axis = sf.axis(i)
        theory = numarray.zeros([len(axis)],numarray.Float)
        const = 100000
        for j in range(0,len(axis)):
            if abs(axis[j]) < R:
                theory[j] = Q/(2.0*R**3)*(3*R**2 - axis[j]**2)
            else:
                theory[j] = Q/abs(axis[j])
        pylab.subplot(2,2,i+1)
        pylab.plot(axis,theory,'r')

if __name__ == "__main__":
    import arrayio

    pylab.figure(1)
    rho = readscalarfield("rho")
    phi = readscalarfield("phi")
#     display_axes(phi)
#     display_theory(phi)
#     pylab.figure(2)
#     display_axes(rho)
    display323232(rho)

    pylab.figure(2)
    display323232(phi)
    
    pylab.show()
