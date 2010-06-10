#!/usr/bin/env python

class Computational_grid:
    def __init__(self, nx, ny, nz, boundary_conditions):
        self.space_charge_bcs = {"3d open":1,
                                 "trans open, long periodic":2,
                                 "trans finite, long open round":3,
                                 "trans finite, long periodic round":4,
                                 "trans finite, long open rect":5,
                                 "trans finite, long periodic rect":6}
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.bc_num = self.space_charge_bcs[boundary_conditions]

    def get_nx(self):
        return self.nx
    def get_ny(self):
        return self.ny
    def get_nz(self):
        return self.nz
    def get_bc_num(self):
        return self.bc_num

