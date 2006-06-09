#!/usr/bin/env python


import sys,os
sys.path.append("..")

import pylab
import diagnostics
import time

#     b  : blue
#     g  : green
#     r  : red
#     c  : cyan
#     m  : magenta
#     y  : yellow
#     k  : black
#     w  : white
colors = ['b','g','r','c','m','y','k']
tuples = []


for i in range(1,len(sys.argv)):
    success = 0
    while not success:
        try:
            d = diagnostics.Diagnostics(sys.argv[i])
            if len(d.num_part) > 0:
                success = 1
        except:
            print "failed to open",sys.argv[i]," will retry in 5 seconds"
            print "queue:"
            os.system("jqstat")
            success = 0
            time.sleep(5)
    tuples.append((d,colors[(i-1) % len(colors)],sys.argv[i]))

#pylab.ion()
pylab.subplot(2,3,1)
for (d,color,name) in tuples:
    pylab.plot(d.s/474.2,d.mean[:,diagnostics.x],color)
pylab.title('x mean')

pylab.subplot(2,3,4)
for (d,color,name) in tuples:
    pylab.plot(d.s/474.2,d.std[:,diagnostics.x],color)
pylab.title('x width')


pylab.subplot(2,3,2)
for (d,color,name) in tuples:
    pylab.plot(d.s/474.2,d.mean[:,diagnostics.y],color)
pylab.title('y mean')

pylab.subplot(2,3,5)
for (d,color,name) in tuples:
    pylab.plot(d.s/474.2,d.std[:,diagnostics.y],color)
pylab.title('y width')


pylab.subplot(2,3,3)
for (d,color,name) in tuples:
    pylab.plot(d.s/474.2,d.num_part/max(d.num_part),color)
limits = pylab.axis()
limits[2] = 0.0
limits[3] = 1.1
pylab.axis(limits)
pylab.title("normalized current")

pylab.subplot(2,3,6)
for (d,color,name) in tuples:
    pylab.plot([0,1],[0,0],color,label=name,linewidth=4.0)
limits[2] = 0.0
limits[3] = 1.0
pylab.axis(limits)
pylab.legend()

pylab.show()


