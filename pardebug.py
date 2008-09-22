#!/usr/bin/env python

from mpi4py import MPI

debug = True
if debug:
    debugfile = open("debug-%02d" % MPI.rank, "w")
else:
    debugfile = None

def pardebug(message):
    if debugfile:
        debugfile.write(message)
        debugfile.flush()
