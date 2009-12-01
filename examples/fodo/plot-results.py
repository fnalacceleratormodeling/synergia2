#!/usr/bin/env python

from matplotlib import pyplot
from diagnostics_file import Diagnostics_file

last = Diagnostics_file("fodo.h5")
zero_current = Diagnostics_file("fodo-np0.h5")
pyplot.plot(last.s,last.std[:,0])
pyplot.plot(zero_current.s,zero_current.std[:,0])
pyplot.show()
