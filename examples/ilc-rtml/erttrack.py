#!/usr/bin/env bwpython

import local_paths
import beamline
import gourmet
import numpy
import physics_constants
import bunch
#import diagnostics
try:
    import pylab
except:
    print "Warning: could not import pylab"
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
#import field
import math

import sys
import memory

import job_manager

import UberPkgpy #SpaceChargePkgpy

import syn2_diagnostics

import fish_fftw as fish
import fish_gauss as fish2
import macro_bunch

import error_eater
import options

import tracker

from mpi4py import MPI

import glob

import loadfile

t2i = {}
t2i_index = -1
def type2int(tname):
    global t2i
    global t2i_index
    if not t2i.has_key(tname):
        t2i_index += 1
        t2i[tname] = t2i_index
    return t2i[tname]
    
def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return int(retval)

def read_em(index,units):
    x = loadfile.loadfile_transpose("code_parts_%d.dat" % index)
    shape = numpy.shape(x)
    for i in range(0,6):
        for j in range(0,shape[1]):
            x[i,j] *= units[i]
    return x

def compare(bunch,index,units,typenum):
    leo = read_em(index,units)
    syn = bunch.get_local_particles()[0:6,0:10]
    passes = True
    sum = 0.0
    for i in range(0,6):
        for j in range(0,10):
            if leo[i,j] == 0.0:
                norm = 1.0
            else:
                norm = leo[i,j]
            norm_diff = (leo[i,j] - syn[i,j])/norm
            sum += norm_diff
            if (abs(norm_diff) > 1.0e-5):
                print i,j,norm_diff,leo[i,j],syn[i,j]
                passes = False
    print "average diff =", sum/60.0,typenum
    return passes

def read_parts(bunch,index,units):
    x = read_em(index,units)
    bunch.get_local_particles()[0:6,0:10] = x

mem0 = 0.0

if ( __name__ == '__main__'):
    t0 = time.time()

#
# Beam parameters at entrance:
# TWISS_EGETAWAY : BETA0, BETX = 1.6479, ALFX =  0.4982, &
##                          BETY = 8.8630, ALFY = -2.3771
#
##  BEAM_EGETAWAY  :  BEAM, ENERGY = 5, NPART = BunchCharge, &
##                    SIGE = 1.5E-3, SIGT = 9.E-03, &
##                    EX = 8E-6 * EMASS/5, EY = 20E-9 * EMASS/5
# where
# EMASS = 0.0005110034
# BunchCharge =  2e10
# where average Current is N/L_scale*q*c
# sigE de/e
# sigt = bunch length (sigma) in meters
# Ex/y = geometric emittance x/y



    myopts = options.Options("ertml")
    myopts.add("current",40.,"current",float)
    # transverse = 1 for longitudinally uniform beam
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("open_all",1,"open transverse boundary conditions", int)
    myopts.add("sigmaz_m",9.E-03,"sigma z in meters",float)
    myopts.add("dpop",1.50e-3,"(delta p)/p",float)
    myopts.add("showplot",0,"show plot",int)
    myopts.add("kicksperline",1000,"kicks per line",int)
    myopts.add("xinit",0.0,"x initial",float)
    myopts.add("xpinit",0.0,"x prime initial",float)
    myopts.add("yinit",0.0,"y initial",float)
    myopts.add("ypinit",0.0,"y prime initial",float)
    myopts.add("dpinit",0.0,"dp initial",float)
    myopts.add("scgrid",[17,17,17],"Space charge grid",int)
    myopts.add("part_per_cell",1,"particlels per cell",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[2,7],"fraction of particles to track (numer,denom)",int)
    myopts.add("xsif","ertml_filecalls.xsif",str)
    myopts.add("exact_propagate",0,int)
    
    myopts.add_suboptions(job_manager.opts)
    myopts.parse_argv(sys.argv)

    job_mgr = job_manager.Job_manager(sys.argv,myopts,
                                      glob.glob("*.xsif"))
    
    f = open("command","w")
    for arg in sys.argv:
        f.write("%s "%arg)
    f.write("\n")
    f.close()

    scaling_frequency = 1.30e9

    pipe_radius = 0.01
    griddim = myopts.get("scgrid")
    part_per_cell = myopts.get("part_per_cell")
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,MPI.COMM_WORLD.Get_size())
    num_particles = 10

##    num_particles = adjust_particles(10,MPI.COMM_WORLD.Get_size())

    ee = error_eater.Error_eater()
    ee.start()

# Need to set positrons!

    g = gourmet.Gourmet(myopts.get("xsif"),"ERTML",5.0-physics_constants.PH_NORM_me,
                        scaling_frequency,myopts.get("maporder"),particle='positron')
    #~ g.insert_element_space_charge_markers(1)
    units = g.get_u(g.get_initial_energy())
###    print units
###    sys.exit(1)

    EMASS = 0.0005110034
    Ex = 8E-6 * EMASS/5
    Ey = 20E-9 * EMASS/5
    rtbetax=1.6479
    rtbetay=8.8630
    alphax =  0.4982
    alphay = -2.3771
    [sigma_x,sigma_xprime,r_x]=matching.match_twiss_emittance(Ex,alphax,rtbetax)   
    [sigma_y,sigma_yprime,r_y]=matching.match_twiss_emittance(Ey,alphay,rtbetay)
    print " twiss matching ", sigma_x, sigma_y, sigma_xprime, sigma_yprime, r_x, r_y

# Is it enough to set the positron mass here?
    bp = beam_parameters.Beam_parameters(physics_constants.PH_NORM_me,
                                         1.0, 5.0, 0.0, scaling_frequency,
                                         myopts.get("transverse"))
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = sigma_x, lam = sigma_xprime * pz)
    bp.y_params(sigma = sigma_y, lam = sigma_yprime * pz)
    sigma_z_meters = myopts.get("sigmaz_m")
    bp.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)
    bp.correlation_coeffs(xpx = r_x, ypy = r_y)


##    bp.x_params(sigma = myopts.get("xinit"), lam = myopts.get("xpinit")*pz)
##    bp.y_params(sigma = myopts.get("yinit"), lam = myopts.get("ypinit")*pz)
##    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
##                     scaling_frequency/math.pi * myopts.get("zfrac")
##    bp.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)
##    bp.correlation_coeffs(xpx = 0.0, ypy = 0.0)

    pgrid = processor_grid.Processor_grid(1)

##    if myopts.get("open_all")==1:
##        cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
##                                                      "3d open")
##    else:
##        cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
##                                                      "trans finite, long periodic round")
    piperad = 0.04
    #field = field.Field(bp, pgrid, cgrid, piperad)
    #
    # Setup old bunch
    #
    old_bunch = bunch.Bunch(myopts.get("current"), bp, num_particles, pgrid)
    old_bunch.generate_particles()

    ####b.write_particles("initial.dat")
     
    line_length = g.orbit_length()
    tau = 0.5*line_length/myopts.get("kicksperline")
    s = 0.0
    ###b.write_fort(s)
    if myopts.get("track"):
        mytracker = tracker.Tracker('/tmp',myopts.get("trackfraction"))


    line_x = None
    line_y = None
    steps = 0
    #~ if MPI.COMM_WORLD.Get_rank() == 0 and myopts.get("showplot"):
        #~ pylab.ion()
        #~ pylab.hold(0)
        #~ xpl=[]
        #~ xppl=[]
        #~ ypl=[]
        #~ yppl=[]
        #~ xx=[]
        #~ yy=[]
        #~ rr=[]
        #~ #
        #~ # Set the firt particle for plotting purposes
        #~ #
        #~ b.particles()[0,0]=myopts.get("xinit")*units[0]
        #~ b.particles()[1,0]=myopts.get("xpinit")*units[1]
        #~ b.particles()[2,0]=myopts.get("yinit")*units[2]
        #~ b.particles()[3,0]=myopts.get("ypinit")*units[3]
        #~ b.particles()[4,0]=0.0
        #~ b.particles()[5,0]=myopts.get("dpinit")*units[5]
#
# Now initialize new bunch
#

    b = macro_bunch.Macro_bunch(physics_constants.PH_NORM_me,1)
    b.init_from_bunch(old_bunch)
    
    read_parts(b,-1,units)
    # and diagnostics
    d = syn2_diagnostics.Diagnostics(units)
    d.add(0,b)

    print "init passes:",compare(b,-1,units,type2int("none"))
    for step in range(0,g.get_num_elements()):
        #~ print "step",step,"= map"
        if myopts.get("exact_propagate") == "1":
            print "exact propagate!"
            g.propagate_element(step,b)
        else:
            print "map propagate 8-(",myopts.get("exact_propagate"),type(myopts.get("exact_propagate"))
            g.get_element_fast_mapping(step).apply(b.get_local_particles(), b.get_num_particles_local())
        s += g.get_element_length(step)        
        print step,g.get_element(step).Name(),g.get_element(step).Type()
        if (step>0) and (step<4613):
            print compare(b,step-1,units,type2int(g.get_element(step).Type()))
            read_parts(b,step-1,units)
        d.add(s,b)
    print "elapsed time =",time.time() - t0
    d.write_hdf5("ertml_fish")

    for key in t2i.keys():
        print key,t2i[key]

