#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import MLab
import physics_constants
import bunch
import diagnostics
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
import field
import math
import basic_toolkit
import beamline

import sys
import memory
from math import pi
from math import sqrt
import octapy


import UberPkgpy #SpaceChargePkgpy

import chef_propagate
import error_eater
import options
import job_manager
import tracker


from mpi4py import MPI
try:
    import pylab
    import density_plot
except:
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "pylab not available"

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return retval

def correct_for_dispersion(b,the_map):
    for i in range(0,b.num_particles_local()):
         b.particles()[0,i] += b.particles()[5,i]*the_map[0,5]

class Line:
    def __init__(self, mad_file, line_name, kinetic_energy, scaling_frequency,
                 map_order, kicks, opts, particle='proton'):
        self.name = line_name
        self.scaling_frequency = scaling_frequency
        self.kicks = kicks
        self.opts = opts
        self.particle = particle
        self.g = gourmet.Gourmet(mad_file,line_name, kinetic_energy,
                                 scaling_frequency, map_order, particle)
        self.g.insert_space_charge_markers(kicks)
	self.g.check()
	self.g.get_single_chef_mapping()
        self.line_length = self.g.orbit_length()
        self.tau = 0.5*self.line_length/self.kicks

    def get_mass(self):
        if self.particle == 'positron':
            return physics_constants.PH_NORM_me
        else:
            return physics_constants.PH_NORM_mp
        
    def propagate(self, s, bunch, sc_params):
        if MPI.COMM_WORLD.Get_rank() == 0:
            print "%s" % self.name,
        mean = None
        std = None
        sys.stdout.flush()
        first_action = 1
        for action in self.g.get_actions():
            if action.is_mapping():
                action.get_data().apply(bunch.particles(),
                                   bunch.num_particles_local())
                s += action.get_length()
            elif action.is_synergia_action():
                if action.get_synergia_action() == "space charge endpoint":
                    if not first_action:
                        t0 = time.time()
                        print "write_fort...",
                        sys.stdout.flush()
                        (mean,std) = b.write_fort(s)
                        t1 = time.time()
                        print " took",t1-t0,"s"
                elif action.get_synergia_action() == "space charge kick":
                    UberPkgpy.Apply_SpaceCharge_external(
                        bunch.get_beambunch(),
                        sc_params.pgrid.get_pgrid2d(),
                        sc_params.field.get_fieldquant(),
                        sc_params.field.get_compdom(),
                        sc_params.field.get_period_length(),
                        sc_params.cgrid.get_bc_num(),
                        sc_params.field.get_pipe_radius(),
                        self.tau, 0, self.scaling_frequency,0)
                else:
                    print "Whachou talkin' about, Willis? '%s'" % \
                          action.get_synergia_action()
            else:
                print "action",action.get_type(),"unknown"
            first_action = 0
        return (s,mean,std)

def get_beam_parameters_envelope_match(line, opts):
    bp = beam_parameters.Beam_parameters(line.get_mass(),
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         opts.get("transverse"))
    [sigma_x,sigma_xprime,r_x,\
     sigma_y,sigma_yprime,r_y] = matching.envelope_match(
        opts.get("emittance"),opts.get("current"),line.g)
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    mismatchx = 1.0 + opts.get("mismatchfracx")
    mismatchy = 1.0 + opts.get("mismatchfracy")
    bp.x_params(sigma = sigma_x*mismatchx, lam = sigma_xprime/mismatchx * pz)
    bp.y_params(sigma = sigma_y*mismatchy, lam = sigma_yprime/mismatchy * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
                     scaling_frequency/math.pi * opts.get("zfrac")
    bp.z_params(sigma = sigma_z_meters, lam = opts.get("dpop")* pz,
                num_peaks = 5)
    bp.correlation_coeffs(xpx = -r_x, ypy = -r_y)

    return bp

def get_beam_parameters_lattice_fns(line, opts):
    bp = beam_parameters.Beam_parameters(line.get_mass(),
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         opts.get("transverse"))
    (alpha_x, alpha_y, beta_x, beta_y) = matching.get_alpha_beta(line.g)
    (sigma_x, sigma_xprime,r_x) = matching.match_twiss_emittance(
        opts.get("emittance"),alpha_x,beta_x)
    (sigma_y, sigma_yprime,r_y) = matching.match_twiss_emittance(
        opts.get("emittance"),alpha_y,beta_y)
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = sigma_x, lam = sigma_xprime * pz)
    bp.y_params(sigma = sigma_y, lam = sigma_yprime * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
                     scaling_frequency/math.pi * opts.get("zfrac")
    bp.z_params(sigma = sigma_z_meters, lam = opts.get("dpop")* pz,
                num_peaks = 1)
    bp.correlation_coeffs(xpx = r_x, ypy = r_y)

    return bp

class Sc_params:
    def __init__(self, griddim, bp, proccol):
        self.pgrid = processor_grid.Processor_grid(proccol)
        self.cgrid = computational_grid.Computational_grid(griddim[0],
                                                      griddim[1],
                                                      griddim[2],
                                                      "trans finite, long periodic round")
        piperad = 0.04
        self.field = field.Field(bp, self.pgrid, self.cgrid, piperad)



if ( __name__ == '__main__'):
    t0 = time.time()

    myopts = options.Options("syn2booster")
    myopts.add("current",0.035,"current",float)
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",3.0e-6,"emittance",float)
    myopts.add("zfrac",0.1,"z width as fraction of bucket",float)
    myopts.add("dpop",5.0e-4,"(delta p)/p",float)
    myopts.add("showplot",0,"show plot",int)
    myopts.add("plotperiod",10,"how often to plot")
    myopts.add("kicksperturn",2000,"kicks per cell",int)
    myopts.add("plotperiod",4,"update plot every plotperiod steps",int)
    myopts.add("turns",12,"number of booster revolutions",int)
    myopts.add("injturns",12,"myoptsnumber of injection turns",int)
    myopts.add("scgrid",[33,33,33],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[2,7],"fraction of particles to track (numer,denom)",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)
    myopts.add("mismatchfracx",0.0,"fractional horizontal mismatch",float)
    myopts.add("mismatchfracy",0.0,"fractional vertical mismatch",float)
    myopts.add("proccol",2,"number of columns in processor grid (y direction)",int)

    myopts.add_suboptions(job_manager.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = job_manager.Job_manager(sys.argv,myopts,
                                      ["ocs6-7.lat",
                                       "envelope_match.cache"])
### start save
    scaling_frequency = 37.7e6
### end save
    part_per_cell = myopts.get("partpercell")
    pipe_radius = 0.04
    griddim = myopts.get("scgrid")
    proccol = myopts.get("proccol")
    if MPI.COMM_WORLD.Get_size() == 1:
        proccol = 1
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,MPI.COMM_WORLD.Get_size())

    ee = error_eater.Error_eater()
#    ee.start()
    cell_line = range(0,25)
    order = myopts.get("maporder")
    kicks = myopts.get("kicksperturn")

    kinetic_energy = sqrt(5.066**2 - physics_constants.PH_NORM_me**2) - \
                     physics_constants.PH_NORM_me
    dr = Line("ocs6-7.lat","ring",kinetic_energy,scaling_frequency,order,
              kicks,myopts,'positron')

    bp = get_beam_parameters_lattice_fns(dr,myopts)

    sc_params = Sc_params(griddim,bp,proccol)
    b = bunch.Bunch(myopts.get("current"), bp, num_particles, sc_params.pgrid)
    b.generate_particles()
    linear_map = dr.g.get_single_linear_map()
    correct_for_dispersion(b,linear_map)

    #b.write_particles("initial.dat")

    last_inj_turn = myopts.get("injturns")
    last_turn = myopts.get("turns")

    s = 0.0
    mean = Numeric.zeros(6,'d')
    std = None
#    b.write_fort(s)

    if myopts.get("track"):
        mytracker = tracker.Tracker('/scratch',myopts.get("trackfraction"))
    time_file = open("timing.dat","w")
    t1 = time.time()
    for turn in range(1,last_turn+1):
        if MPI.COMM_WORLD.Get_rank() == 0:
            time_file.write("%g %g\n" % (s,time.time()-t1))
            t1 = time.time()
            time_file.flush()
        if MPI.COMM_WORLD.Get_rank()==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        print "jfa: propagate start"
        (s,mean,std) = dr.propagate(s,b, sc_params)
        print "jfa: propagate end"
        (mean,std) = b.write_fort(s)
        if MPI.COMM_WORLD.Get_rank()==0:
            print        
        if turn % myopts.get("saveperiod") == 0:
            b.write_particles("turn_%04d.dat" % turn)
    if myopts.get("track"):
        mytracker.close()
        mytracker.show_statistics()
    MPI.WORLD.Barrier()
    print "elapsed time =",time.time() - t0
    time_file.close()

