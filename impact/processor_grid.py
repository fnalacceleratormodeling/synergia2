#!/usr/bin/env python

from mpi4py import MPI
from UberPkgpy import *

class Processor_grid:
    def __init__(self,columns):
        self.pgrid2d = Pgrid2d()
        self.npcol = columns
        self.nprow = MPI.size/self.npcol
        mpi_comm_world = fmpi_comm_world_external()
        construct_Pgrid2d_external(self.pgrid2d, fmpi_comm_world_external(),
                                   self.nprow, self.npcol)

    def get_pgrid2d(self):
        return self.pgrid2d
    def get_nrow(self):
        return self.nprow
    def get_ncol(self):
        return self.npcol
