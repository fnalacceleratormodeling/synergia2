#!/usr/bin/env python

from mpi4py import MPI

debugfile = None

def pardebug(message):
    global debugfile
    if not debugfile:
        debugfile = open("debug-%02d" % MPI.COMM_WORLD.Get_rank(), "w")
    debugfile.write(message)
    debugfile.flush()
