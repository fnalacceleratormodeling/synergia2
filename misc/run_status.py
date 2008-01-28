#!/usr/bin/env bwpython

import diagnostics
import pylab
import time

pylab.ion()

while 1:
    found = 0
    try:
        d = diagnostics.Diagnostics("run1")
        found = 1
    except:
        print "failed to open diagnostics, will try again in 30 seconds"
    if found:
        line_x, = pylab.plot(d.s/474.2,d.std[0:96:d.num,diagnostics.x],'r-o',label='x')
        line_y, = pylab.plot(d.s/474.2,d.std[0:96:d.num,diagnostics.y],'b-x',label='y')
        pylab.draw()
    time.sleep(30)
