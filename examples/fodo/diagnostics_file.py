#!/usr/bin/env python

import tables

class Diagnostics_file:
    def __init__(self,filename):
        d = tables.openFile(filename)
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
