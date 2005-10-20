#!/usr/bin/env python

from mpi4py import MPI
from Pgrid2dPkgpy import *

# Oi. We need a more general solution for this
MPI.MPI_COMM_WORLD=91

class Processor_grid:
    def __init__(self):
        self.pgrid2d = Pgrid2d()
        # fixme
        self.nprow = 1
        self.npcol = 1
        construct_Pgrid2d_external(self.pgrid2d, MPI.MPI_COMM_WORLD,
                                   self.nprow, self.npcol)

    def get_pgrid2d(self):
        return self.pgrid2d
    def get_nrow(self):
        return self.nprow
    def get_ncol(self):
        return self.npcol
