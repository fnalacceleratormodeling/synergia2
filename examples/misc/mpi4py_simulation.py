#!/home/dechow/Python/bin/bwpython

from mpi4py import MPI
mpi_rank = MPI.rank
mpi_size = MPI.size
print "I am processor %d of %d" % (mpi_rank, mpi_size)

import local_paths
import array
import sys
import Numeric

# FORTHON PACKAGE IMPORTS
from Pgrid2dPkgpy      import *
from BeamBunchPkgpy    import *
from ExternalBLEPkgpy  import *
from DistributionPkgpy import *
from OutputPkgpy       import *

def Create_BeamBunch( incurr, inkin, inmass, incharge, innp, phasini ):
    # THIS IS NOT A REAL CONSTRUCTOR CALL; KEEP PYTHON TYPE SYSTEM HAPPY
    bunch = BeamBunch()
    # NOW WE MAKE A BeamBunch
    construct_BeamBunch_external(bunch, incurr, inkin, inmass, incharge, innp, phasini)
    return bunch

def Create_Pgrid2d(  MPI_COMM_WORLD, nprow, npcol ):
    # AGAIN; NOT A REAL PGrid2d
    grid = Pgrid2d( )
    construct_Pgrid2d_external( grid, MPI_COMM_WORLD, nprow, npcol )
    return grid

def Create_BeamLineElement( bnseg, bmpstp, bitype, blength, frequency, MapOrder ):
    # FOR NOW--9 Jun 05--bury this initialization of mad file input here
    init_mad_external()
    element = ExternalBLE()
    construct_ExternalBLE_external( element, bnseg, bmpstp, bitype, blength, frequency, MapOrder )
    # ARE THE NEXT TWO LINES PART OF THE CONSTRUCTION PROCESS?
    tmpextble = Numeric.array([ 0.0, 1.0, 0., 0.040000, 1.0, 0., 0., 0., 0. ] ) 
    setparam2_ExternalBLE_external( element, tmpextble )
    print "jfa: external_ble.mapstp =", element.Mapstp
    return element

def Initial_Distribution( bunch, nparam, distparam, grid, flagalloc ):
    Gauss_Covariance_Dist_external( bunch, nparam, distparam, grid, flagalloc )
    return bunch

def Apply_Map(a_beam_bunch, a_map_1, a_map_2 = None):
    if ( a_map_2 ):
	MapOrder = 2
    else:
	MapOrder = 1
    m = Numeric.zeros( (7,7), 'd' )
    m[0:6,0:6] = a_map_1
    m[6,6] = 1.0
    a_beam_bunch.Pts1 = Numeric.matrixmultiply( m, a_beam_bunch.Pts1 )

def Run_Simulation( bunch, element, xm, z, tau, number_of_turns, number_of_segments ):
    # CREATE A DATA FILE FOR VERIFICATION IN OCTAVE
    mapfile = open("map.dat","w")

    #diagnostic3_Output_external( z, bunch, frequency )
    # note the difference--no external proxy to call through
    diagnostic3_Output( z, bunch, frequency )
    # go for number_of_turns
    for turns in range( 0, number_of_turns ):
	# go for number_of_segments
	for segments in range( 0, number_of_segments ):
	    maplinear_ExternalBLE_external((z - turns * element.Length), tau, xm, element, bunch.refptcl, bunch.Charge, bunch.Mass )
	    # FILL A DATA FILE FOR VERIFICATION IN OCTAVE
	    for row in range(0,6):
		for col in range(0,6):
		    mapfile.write( "%g " % xm[row,col] )
		mapfile.write( "\n" )

	    Apply_Map( bunch, xm )
	    z = z + tau
	    maplinear_ExternalBLE_external((z - turns * element.Length), tau, xm, element, bunch.refptcl, bunch.Charge, bunch.Mass )

	    # FILL A DATA FILE FOR VERIFICATION IN OCTAVE
	    for row in range(0,6):
		for col in range(0,6):
		    mapfile.write( "%g " % xm[row,col] )
		mapfile.write( "\n" )

	    Apply_Map( bunch, xm )
	    z = z + tau
	    #diagnostic3_Output_external( z, bunch, frequency )
	    # note the difference--no external proxy to call through
	    diagnostic3_Output( z, bunch, frequency )

    # CLOSE THE DATA FILE FOR VERIFICATION IN OCTAVE
    mapfile.close( )

if ( __name__ == '__main__'):
    # simulation parameters
    # Beam bunch paramenters
    incurr = 0                     # obcurr in Input.f90
    inkin = 4e+8                   # obkenergy
    inmass = 9.38272e+8            # obmass
    incharge = 1.0                 # obcharge
    phasini = 0.0                  # ophsini
    innp = 1000                    # onp

    # Grid parameters
    MPI_COMM_WORLD = 91
    nprow = 1                      # onprow
    npcol = 2                      # onpcol
    time =  0.0

    # Distribution parameters
    nparam = 30
    flagalloc = 0
#    distparam = array.array( 'd', [ 0.000281124, -2.25491e-20, 4.9945e-07, 0, 0, 0.000991052, 0,
#				    0, -1.55282e-20, 1.41675e-07, 0, 0, 0, 0,
#				    0.345566, -6.62671e-07, 0, 0, 0, 0, 9.30938e-08, 0, 0, 0, 0, 0, 0, 1, 666, 666] )
    distparam = Numeric.array( [ 0.000281124, 0, 4.9945e-07, 0, 0, 0.000991052, 0,
				 0, 0, 1.41675e-07, 0, 0, 0, 0,
				 0.345566, -6.62671e-07, 0, 0, 0, 0, 9.30938e-08, 0, 0, 0, 0, 0, 0, 1, 666, 666] ,'d' )

    # External Beam Line Element Parameters
    bnseg = 100
    bmpstp = 1
    bitype = 91
    blength = 474.203
    MapOrder = 1
    frequency = 2e+08
    
    # Map parameters
    xm = Numeric.zeros( (6,6), 'd' )
    z = 0.0
    tau = 0.5*blength/bnseg
    number_of_turns = 10

    # initialize the necessary simulation components
    the_bunch = Create_BeamBunch( incurr, inkin, inmass, incharge, innp, phasini )
    the_grid = Create_Pgrid2d(  MPI_COMM_WORLD, nprow, npcol )
    the_bunch = Initial_Distribution(  the_bunch, nparam, distparam, the_grid, flagalloc )
    the_element = Create_BeamLineElement( bnseg, bmpstp, bitype, blength, frequency, MapOrder )

    # run the simulation
    Run_Simulation( the_bunch, the_element, xm, z, tau, number_of_turns, bnseg )

    print "Numeric.shape( the_bunch.Pts1 ): ",Numeric.shape( the_bunch.Pts1 )

    x = MPI.WORLD[ 0 ].Gather( the_bunch.Pts1[0,:] )
    y = MPI.WORLD[ 0 ].Gather( the_bunch.Pts1[1,:] )

    if mpi_rank == 0:
	print "shape x: ", Numeric.shape( x )
	print x
	print "shape y: ", Numeric.shape( y )
	print y
	#file = open("pympi-x1.txt", 'w')
	#for i in range(0,1000):
        #    file.write("%g "% x[i] )
	#    file.write("\n")

#    from pylab import *
#    if mpi_rank == 0:
#	scatter( x, y )
#	show()

   # MPI.Finalize()

