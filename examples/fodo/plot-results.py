#!/usr/bin/env python

import tables
from matplotlib import pyplot

class Diagnostics:
    def __init__(self,filename):
        d = tables.openFile(filename)
#        print dir(d.root)
        self.s = d.root.s.read()
        self.emitx = d.root.emitx.read()
        self.emitxy = d.root.emitxy.read()
        self.emitxyz = d.root.emitxyz.read()
        self.emity = d.root.emity.read()
        self.emitz = d.root.emitz.read()
        self.mean = d.root.mean.read()
        self.n = d.root.n.read()
        self.s = d.root.s.read()
        self.std = d.root.std.read()
        self.units = d.root.units.read()
        d.close()

last = Diagnostics("fodo.h5")
zero_current = Diagnostics("fodo-np0.h5")
pyplot.plot(last.s,last.std[:,0])
pyplot.plot(zero_current.s,zero_current.std[:,0])
pyplot.show()
