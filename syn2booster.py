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


import UberPkgpy #SpaceChargePkgpy

import apply_map
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
    if MPI.rank == 0:
        print "pylab not available"

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return retval

mem0 = 0.0

plot_index = 0
def plot_long(b):
    if MPI.rank ==0:
        z = b.particles()[4,:]
        zprime = b.particles()[5,:]

        d = diagnostics.Diagnostics()
        pylab.hold(0)
        pylab.ioff()

        sample_period = 1
        pylab.subplot(2,2,1)
        pylab.hold(0)
        pylab.plot(d.s[1:len(d.mean):sample_period]/474.2,d.mean[1:len(d.mean):sample_period,diagnostics.x])
        pylab.title("horizontal position")
        pylab.subplot(2,2,2)
        pylab.hold(0)
        pylab.plot(d.s[1:len(d.mean):sample_period]/474.2,
                   d.std[1:len(d.mean):sample_period,diagnostics.x])
        pylab.title("horizontal width")

        pylab.subplot(2,2,4)
        pylab.hold(0)
        pylab.plot(d.s[1:len(d.mean):sample_period]/474.2,
                   d.std[1:len(d.mean):sample_period,diagnostics.xprime])
        pylab.title("horizontal angular width")
#         pylab.plot(d.s/474.2,d.num_part,'g')
#         pylab.title("number of particles")

        pylab.subplot(2,2,3)
        density_plot.density_plot(z,zprime,70)
        pylab.title("phase space")
        limits = pylab.axis()
        limits[0] = -pi
        limits[1] = pi
        pylab.axis(limits)
        pylab.title("long ps")

        pylab.ion()
        pylab.draw()
        global plot_index
#        pylab.savefig("plot%03d.png" % plot_index)
        plot_index += 1
    
def correct_for_dispersion(b,the_map):
    for i in range(0,b.num_particles_local()):
         b.particles()[0,i] += b.particles()[5,i]*the_map[0,5]

rfcells = [14,15,16,17,19,21,22,23,24]
class Line:
    def __init__(self, mad_file, line_name, kinetic_energy, scaling_frequency,
                 map_order, kicks, opts):
        global rfcells
        self.name = line_name
        self.scaling_frequency = scaling_frequency
        self.kicks = kicks
        self.opts = opts
        self.rfnames = ["bcel%02d" % cell for cell in rfcells]
        self.g = gourmet.Gourmet(mad_file,line_name, kinetic_energy,
                                 scaling_frequency, map_order)
        if line_name in self.rfnames:
            self.insert_rf()
        self.g.insert_space_charge_markers(kicks)
        self.line_length = self.g.orbit_length()
        self.tau = 0.5*self.line_length/self.kicks
        
    def insert_rf(self):
        self.g.iterator.reset()
        element = self.g.iterator.next()
        s = 0.0
        while element.Name() != "LONGA:":
            s += element.OrbitLength(self.g.get_initial_particle())
            element = self.g.iterator.next()
        rf_total_length = element.OrbitLength(self.g.get_initial_particle())
        cavity_length = 2.56
        drift_length = (rf_total_length - 2.0*cavity_length)/2.0
        s_rf1 = s + drift_length + 0.5*cavity_length
        s_rf2 = s_rf1 + cavity_length
        phase1 = self.opts.get("rfphase") - math.pi/2.0
        phase2 = self.opts.get("rfphase") + math.pi/2.0
        voltage = self.opts.get("rfvoltage")
        freq = self.scaling_frequency/(2.0*math.pi)
        quality = 0
        impedance = 0
        elements = (beamline.thinrfcavity("synergia action:rfcavity1",
                                          freq,voltage,phase1,quality,impedance),
                    beamline.thinrfcavity("synergia action:rfcavity2",
                                          freq,voltage,phase2,quality,impedance))
        positions = (s_rf1,s_rf2)
        self.g.insert_elements(elements,positions)

    def set_rf_params(self,voltage_ev,phase1,phase2):
        if self.name in self.rfnames:
            self.g.iterator.reset()
            element = self.g.iterator.next()
            while element:
                if element.Type() == 'thinrfcavity':
                    if element.Name() == "synergia action:rfcavity1":
                        phase = phase1
                    elif element.Name() == "synergia action:rfcavity2":
                        phase = phase2
                    else:
                        print "Unknown rf cavity",element.Name()
                    strength = 1.0e-9 * voltage_ev
                    element.setStrength(strength)
                    element.setPhi(phase)
                element = self.g.iterator.next()
            self.g.generate_actions()

    def propagate(self, s, bunch, sc_params):
        if MPI.rank == 0:
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
                        pass
#                         (mean,std) = b.write_fort(s)
                elif action.get_synergia_action() == "space charge kick":
                    UberPkgpy.Apply_SpaceCharge_external(
                        bunch.get_beambunch(),
                        sc_params.pgrid.get_pgrid2d(),
                        sc_params.field.get_fieldquant(),
                        sc_params.field.get_compdom(),
                        sc_params.field.get_period_length(),
                        sc_params.cgrid.get_bc_num(),
                        sc_params.field.get_pipe_radius(),
                        self.tau, 0, self.scaling_frequency)
                elif action.get_synergia_action() == "rfcavity1" or \
                     action.get_synergia_action() == "rfcavity2":
                    element = action.get_data()
                    u_in = self.g.get_u(action.get_initial_energy())
                    u_out = self.g.get_u(action.get_final_energy())
                    chef_propagate.chef_propagate(
                        bunch.particles(), bunch.num_particles_local(),
                        element, action.get_initial_energy(),
                        u_in, u_out)
                else:
                    print "Whachou talkin' about, Willis? '%s'" % \
                          action.get_synergia_action()
            else:
                print "action",action.get_type(),"unknown"
            first_action = 0
        return (s,mean,std)

def get_beam_parameters_orig(line, opts):
    bp = beam_parameters.Beam_parameters(physics_constants.PH_NORM_mp,
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         opts.get("transverse"))
    [sigma_x,sigma_xprime,r_x,\
     sigma_y,sigma_yprime,r_y] = matching.envelope_match(
        opts.get("emittance"),opts.get("current"),line.g)
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = sigma_x, lam = sigma_xprime * pz)
    bp.y_params(sigma = sigma_y, lam = sigma_yprime * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
                     scaling_frequency/math.pi * opts.get("zfrac")
    bp.z_params(sigma = sigma_z_meters, lam = opts.get("dpop")* pz,
                num_peaks = 5)
    bp.correlation_coeffs(xpx = -r_x, ypy = -r_y)

    return bp

def get_beam_parameters(line, opts):
    bp = beam_parameters.Beam_parameters(physics_constants.PH_NORM_mp,
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
    def __init__(self, griddim, bp):
        self.pgrid = processor_grid.Processor_grid(1)
        self.cgrid = computational_grid.Computational_grid(griddim[0],
                                                      griddim[1],
                                                      griddim[2],
                                                      "trans finite, long periodic round")
        piperad = 0.04
        self.field = field.Field(bp, self.pgrid, self.cgrid, piperad)


def update_rf(cell_line,opts,turn,beam_mean):
    global rfcells
    for cell in rfcells:
        if turn < (opts.get("norfturns")+opts.get("rampturns")):
            ramp = (turn-opts.get("norfturns"))/(1.0*opts.get("rampturns"))
            offset = pi/2.0*(1-ramp)
        else:
            offset = 0.0
        voltage = opts.get("rfvoltage")
#        deviation = beam_mean[5]*1.0e6 # convert MeV to eV
#         if voltage > 0.0:
#             rpos_correction = -math.asin(deviation/\
#                                          (2*len(rfcells)*voltage*math.cos(offset)))
#             compaction_correction = -beam_mean[4]*math.pi/180.0
#         else:
#             rpos_correction = 0.0
#             compaction_correction = 0.0
        phase0 = opts.get("rfphase")
        cell_line[cell].set_rf_params(voltage,
                                      phase0 - offset,
                                      phase0 + offset)

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
    myopts.add("kickspercell",4,"kicks per cell",int)
    myopts.add("plotperiod",4,"update plot every plotperiod steps",int)
    myopts.add("turns",12,"number of booster revolutions",int)
    myopts.add("rfvoltage",0.6e6/18.0*0,"rf voltage (MV)",float)
    myopts.add("rfphase",0.0,"rf cavity phase (rad)",float)
    myopts.add("rampturns",200,"number of turns for rf phase ramping",int)
    myopts.add("norfturns",0,"number of turns without any rf BROKEN!",int)
    myopts.add("injturns",12,"myoptsnumber of injection turns",int)
    myopts.add("scgrid",[33,33,33],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[2,7],"fraction of particles to track (numer,denom)",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)

    myopts.add_suboptions(job_manager.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = job_manager.Job_manager(sys.argv,myopts,
                                      ["booster_classic.lat",
                                       "envelope_match.cache"])
### start save
    scaling_frequency = 37.7e6
### end save
    part_per_cell = myopts.get("partpercell")
    pipe_radius = 0.04
    griddim = myopts.get("scgrid")
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,MPI.size)

    ee = error_eater.Error_eater()
    ee.start()
    cell_line = range(0,25)
    order = myopts.get("maporder")
    kicks = myopts.get("kickspercell")
    for cell in range(1,25):
        cell_line[cell] = Line("booster_classic.lat","bcel%02d"%cell ,0.4,
                               scaling_frequency,order,kicks,myopts)
    injcell_line = Line("booster_classic.lat","bcelinj" ,0.4,
                        scaling_frequency,order,kicks,myopts)
    part_kicks = 2
    inja_line = Line("booster_classic.lat","bcel01a" ,0.4,
                     scaling_frequency,order,part_kicks,myopts)
    injb_line = Line("booster_classic.lat","bcel01b" ,0.4,
                     scaling_frequency,order,part_kicks,myopts)

    bp = get_beam_parameters_orig(injcell_line,myopts)

    sc_params = Sc_params(griddim,bp)
    b = bunch.Bunch(myopts.get("current"), bp, num_particles, sc_params.pgrid)
    b.generate_particles()
    linear_map = injcell_line.g.get_single_linear_map()
    correct_for_dispersion(b,linear_map)

    #b.write_particles("initial.dat")

    last_inj_turn = myopts.get("injturns")
    last_turn = myopts.get("turns")

    s = 0.0
    mean = Numeric.zeros(6,'d')
    std = None
#    b.write_fort(s)

    if MPI.rank == 0 and myopts.get("showplot"):
        pylab.ion()
    if myopts.get("track"):
        mytracker = tracker.Tracker('/scratch',myopts.get("trackfraction"))
    time_file = open("timing.dat","w")
    t1 = time.time()
    for turn in range(1,last_inj_turn+1):
        if MPI.rank == 0:
            time_file.write("%g %g\n" % (s,time.time()-t1))
            t1 = time.time()
            time_file.flush()
        update_rf(cell_line,myopts,turn,mean)
        if MPI.rank==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        (s,mean,std) = injb_line.propagate(s, b, sc_params)
        (mean,std) = b.write_fort(s)
        if turn % myopts.get("plotperiod") == 0 and myopts.get("showplot"):
            plot_long(b)
        for cell in range(2,25):
            (s,mean,std) = cell_line[cell].propagate(s,b, sc_params)
            (mean,std) = b.write_fort(s)
            if cell % 24 == 0 and turn % myopts.get("plotperiod") == 0 \
                   and myopts.get("showplot"):
                plot_long(b)
            if cell % 12 == 2 and myopts.get("track"):
                mytracker.add(b,s)
        if turn < last_inj_turn:
            (s,mean,std) = inja_line.propagate(s,b,sc_params)
            binj = bunch.Bunch(myopts.get("current"), bp, num_particles,
                               sc_params.pgrid)
            binj.generate_particles()
            correct_for_dispersion(binj,injcell_line.g.get_single_linear_map())
            b.inject(binj)
        if MPI.rank==0:
            print
        if turn % myopts.get("saveperiod") == 0:
            b.write_particles("turn_inj_%04d.dat" % turn)
    for turn in range(last_inj_turn+1,last_turn+1):
        if MPI.rank == 0:
            time_file.write("%g %g\n" % (s,time.time()-t1))
            t1 = time.time()
            time_file.flush()
        update_rf(cell_line,myopts,turn,mean)
        if MPI.rank==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        for cell in range(1,25):
            (s,mean,std) = cell_line[cell].propagate(s,b, sc_params)
            (mean,std) = b.write_fort(s)
            if cell % 24 == 0 and turn % myopts.get("plotperiod") == 0 \
                   and myopts.get("showplot"):
                plot_long(b)
            if cell % 12 == 2 and myopts.get("track"):
                mytracker.add(b,s)
        if MPI.rank==0:
            print        
        if turn % myopts.get("saveperiod") == 0:
            b.write_particles("turn_%04d.dat" % turn)
    if myopts.get("track"):
        mytracker.close()
        mytracker.show_statistics()
    MPI.WORLD.Barrier()
    print "elapsed time =",time.time() - t0
    time_file.close()

