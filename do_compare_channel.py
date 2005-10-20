#!/usr/bin/env python

import diagnostics
import pylab
import loadfile

def doit():
    d = diagnostics.Diagnostics()

    pylab.plot(d.s,d.std[:,diagnostics.x],'ro')
    d0 = diagnostics.Diagnostics("channel0current")

    pylab.plot(d0.s,d0.std[:,diagnostics.x],'gx')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    e = loadfile.loadfile("/home3/amundson/work/envelope-comparisons/envelope_match_channel_0.5A.dat")
    pylab.plot(e[:,0],e[:,1])

    dold = diagnostics.Diagnostics("cn.06")
    pylab.plot(dold.s,dold.std[:,diagnostics.x],'yo')
    pylab.show()

if ( __name__ == '__main__'):
    doit()
