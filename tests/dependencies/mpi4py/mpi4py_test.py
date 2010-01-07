#!/usr/bin/env python

import os
import sys
from mpi4py_test_options import opts

if __name__ == "__main__":
    print "hello world from mpi4py_test.py"

    print(os.uname())
    sys.stdout.flush()
    
    print(sys.version)
    sys.stdout.flush()
    
#    print(sys.path)
#    sys.stdout.flush()
    
    print(sys.platform)
    sys.stdout.flush()
    
#    print(sys.modules)
#    sys.stdout.flush()
    
    print(sys.argv)
    sys.stdout.flush()
    
    import mpi4py
    import mpi4py.MPI
    
    mysize = mpi4py.MPI.COMM_WORLD.Get_size()
    sys.stdout.write("my size = %d\n" % mysize)
    sys.stdout.flush()
    
    myname = mpi4py.MPI.COMM_WORLD.Get_name()
    sys.stdout.write("my name = %s\n" % myname)
    sys.stdout.flush()
    
    myrank = mpi4py.MPI.COMM_WORLD.Get_rank()
    sys.stdout.write("my rank = %d\n" % myrank)
    sys.stdout.flush()
    
    outnm = "mpiout.%d" % myrank
    
    f = open(outnm, "w")
    
    f.write("my size = %d\n" % mysize)
    f.flush()
    
    f.write("my rank = %d\n" % myrank)
    f.flush()
    
    f.write("my name = %s\n" % myname)
    f.flush()
    
#    dir(mpi4py.MPI.COMM_WORLD)
#    sys.stdout.flush()
    
    sys.stderr.write("before Allgather\n")
    sys.stderr.flush()
    
    gatherres = mpi4py.MPI.COMM_WORLD.allgather(myrank)
    sys.stderr.write("after Allgather\n")
    sys.stderr.flush()
    
    f.write("results of Allgather: %s\n" % repr(gatherres))
    f.flush()
    sys.stdout.write("results of Allgather: %s\n" % repr(gatherres))
    sys.stdout.flush()
    
    
    sumres = mpi4py.MPI.COMM_WORLD.allreduce(myrank, op=mpi4py.MPI.SUM)
    f.write("Results of allreduce(SUM): %s\n" % repr(sumres))
    f.flush()
    sys.stdout.write("Results of Allreduce(SUM): %s\n" % repr(sumres))
    sys.stdout.flush()
