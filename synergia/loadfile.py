#!/usr/bin/env python

import numpy

def loadfile(filename):
    num = len(open(filename ,"r").readlines())
    f = open(filename,"r")
    line = f.readline()
    cols = len(line.split())
    data = numpy.zeros((num,cols),'d')
    index = 0
    while line:
        cols = line.split()
        data[index,:] = map(float,cols)
        line = f.readline()
        index += 1
    f.close()
    return data

def loadfile_transpose(filename):
    num = len(open(filename ,"r").readlines())
    f = open(filename,"r")
    line = f.readline()
    cols = len(line.split())
    data = numpy.zeros((cols,num),'d')
    index = 0
    while line:
        cols = line.split()
        data[:,index] = map(float,cols)
        line = f.readline()
        index += 1
    f.close()
    return data
