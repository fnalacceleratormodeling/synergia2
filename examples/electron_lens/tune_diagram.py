#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from matplotlib import pyplot

def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]

def tune_diagram(max_order,colors='black',linew=1, xbias=0, ybias=0):
    x = numpy.array([0,0],'d')
    y = numpy.array([0,0],'d')
    colors = listify(colors)
    
    # We want to plot
    #   a*x + b*y = p
    # over the range (x,y) = (0,0) ... (1,1)
    # solving for x and y gives
    #   x = (p - b*y)/a
    #   y = (p - a*x)/b
    # the order is given by
    #   ord = a + abs(b)
    for ord in range(max_order,0,-1):
        color = colors[ord%len(colors)]
        for a in range(0,ord+1):
            for sign in range(-1,2,2):
                b = sign*(ord-a)
                if (a == 0):
                #    y = p/b
                #   0 < p < b
                    for p in range(1,b):
                        x[0] = 0
                        y[0] = p/(1.0*b)
                        x[1] = 1
                        y[1] = p/(1.0*b)
                        pyplot.plot(xbias+x,ybias+y,color,linewidth=linew)
                if (b == 0):
                # x = p/a
                # 0 < p < a
                    for p in range(1,a):
                        x[0] = p/(1.0*a)
                        y[0] = 0
                        x[1] = p/(1.0*a)
                        y[1] = 1
                        pyplot.plot(xbias+x,ybias+y,color,linewidth=linew)
                if ((a != 0) and (b !=0)):
                    x[0] = 0
                    x[1] = 1
                    if (b > 0):
                        pmin = 0
                        pmax = a + b
                    else:
                        pmin = b
                        pmax = a
                    for p in range(pmin,pmax+1):
                        y[0] = (p - a*x[0])/(1.0*b)
                        y[1] = (p - a*x[1])/(1.0*b)
                        pyplot.plot(xbias+x,ybias+y,color,linewidth=linew)
                    
if __name__ == "__main__":
    tune_diagram(4,'red')
    pyplot.axis([0,1,0,1])
    pyplot.show()
