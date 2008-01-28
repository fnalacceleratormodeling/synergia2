#!/usr/bin/env bwpython

import diagnostics
import pylab

def doit():
    d = diagnostics.Diagnostics()
    d_nosc = diagnostics.Diagnostics("/home3/amundson/work/Layer-head/transition/cl.01")
    dold = diagnostics.Diagnostics("/home3/amundson/work/Layer-head/transition/bc.03")

    pylab.subplot(211)
    pylab.plot(d.s,d.std[:,diagnostics.x],'r-')
    pylab.plot(dold.s,dold.std[:,diagnostics.x],'g-b')
#    pylab.plot(d_nosc.s,d_nosc.std[:,diagnostics.x],'3')
    pylab.legend(('new', 'old'),'upper center')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    pylab.subplot(212)
    pylab.plot(d.s,d.std[:,diagnostics.xprime],'r-')
    pylab.plot(dold.s,dold.std[:,diagnostics.xprime],'b-')
#    pylab.plot(d_nosc.s,d_nosc.std[:,diagnostics.x],'3')
    pylab.legend(('new', 'old'),'upper center')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<xprime> (m)')

#    pylab.figure()
#    pylab.plot(d.s,d.corr[:,diagnostics.x_xprime])    
#    pylab.plot(d_nosc.s,d_nosc.corr[:,diagnostics.x_xprime])    
#     pylab.ylabel("r_xx' (unitless)")
#     pylab.xlabel('s (m)')
    pylab.show()

if ( __name__ == '__main__'):
    doit()
