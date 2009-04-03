#!/usr/bin/env python

import tables
import pylab
import synergia

#f = tables.openFile("ps2tunes.h5","r")
#tunesx = f.root.tunesx.value.read()
#tunesy = f.root.tunesy.value.read()

#synergia.density_plot(tunesx,tunesy,10)
#f = tables.openFile("ps2tunes_42e11.h5","r")
#f = tables.openFile("ps2tunes_7e11.h5","r")
f = tables.openFile("ps2tunes_7e11_p8000.h5","r")
tunesx = f.root.t1.value.read()
tunesy = f.root.t2.value.read()

#synergia.density_plot(tunesx,tunesy,15)
#text='I=2.69 A\nE=4 GeV\nppb=4.2e+11  '
text='I=4.486 A \nE=4 GeV \nppb=7e+11'
synergia.density_plot_alex(tunesx,tunesy,40,text=text)
pylab.show()
#pylab.savefig('ps2tunes_4.2e11.png',dpi=120) 
#pylab.savefig('ps2tunes_7e11_p8000.png',dpi=150) 