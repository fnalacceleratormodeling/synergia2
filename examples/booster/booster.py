#!/usr/bin/env bwpython

import synergia
import s2_fish

import basic_toolkit
import beamline

import Numeric
import MLab
import time
import math
import sys
from math import pi

from mpi4py import MPI
try:
    import pylab
except:
    if MPI.rank == 0:
        print "pylab not available"
    
def correct_for_dispersion(bunch,the_map):
    for i in range(0,bunch.get_num_particles_local()):
         bunch.get_local_particles()[0,i] += bunch.get_local_particles()[5,i]*the_map[0,5]

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
        self.gourmet = synergia.Gourmet(mad_file,line_name, kinetic_energy,
                                 scaling_frequency, map_order)
        if line_name in self.rfnames:
            self.insert_rf()
        self.gourmet.insert_space_charge_markers(kicks)
        self.line_length = self.gourmet.orbit_length()
        self.tau = 0.5*self.line_length/self.kicks
        self.diag = synergia.Diagnostics(self.gourmet.get_initial_u())
        self.n_sigma = 10
        self.t0 = time.time()
        
    def insert_rf(self):
        s = 0.0
        for element in self.gourmet.beamline:
            if element.Name() == "LONGA":
                break
            s += element.OrbitLength(self.gourmet.get_initial_particle())
        rf_total_length = element.OrbitLength(self.gourmet.get_initial_particle())
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
        self.gourmet.insert_elements(elements,positions)

    def set_rf_params(self,voltage_ev,phase1,phase2):
        if self.name in self.rfnames:
            for element in self.gourmet.beamline:
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
            self.gourmet.generate_actions()

def get_beam_parameters(line, opts):
    beam_parameters = synergia.Beam_parameters(synergia.PH_NORM_mp,
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         opts.get("transverse"))
    [sigma_x,sigma_xprime,r_x,\
     sigma_y,sigma_yprime,r_y] = synergia.envelope_match(
        opts.get("emittance"),opts.get("emittance"),
        opts.get("current"),line.gourmet)
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    mismatchx = 1.0 + opts.get("mismatchfracx")
    mismatchy = 1.0 + opts.get("mismatchfracy")
    beam_parameters.x_params(sigma = sigma_x*mismatchx, 
        lam = sigma_xprime/mismatchx * pz,
        r=-r_x)
    beam_parameters.y_params(sigma = sigma_y*mismatchy, 
        lam = sigma_yprime/mismatchy * pz,
        r=-r_y)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/\
                     scaling_frequency/math.pi * opts.get("zfrac")
    beam_parameters.z_params(sigma = sigma_z_meters, lam = opts.get("dpop")* pz,
                num_peaks = 5)
    return beam_parameters
    
def update_rf(cell_line,opts,turn):
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

    myopts = synergia.Options("booster")
    myopts.add("current",0.035,"current",float)
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",3.0e-6,"emittance",float)
    myopts.add("zfrac",0.1,"z width as fraction of bucket",float)
    myopts.add("dpop",5.0e-4,"(delta p)/p",float)
    myopts.add("kickspercell",4,"kicks per cell",int)
    myopts.add("turns",12,"number of booster revolutions",int)
    myopts.add("rfvoltage",0.6e6/18.0,"rf voltage (MV)",float)
    myopts.add("rfphase",0.0,"rf cavity phase (rad)",float)
    myopts.add("rampturns",200,"number of turns for rf phase ramping",int)
    myopts.add("norfturns",0,"number of turns without any rf BROKEN!",int)
    myopts.add("injturns",12,"myoptsnumber of injection turns",int)
    myopts.add("scgrid",[33,33,33],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[2,7],"fraction of particles to track (numer,denom)",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)
    myopts.add("mismatchfracx",0.0,"fractional horizontal mismatch",float)
    myopts.add("mismatchfracy",0.0,"fractional vertical mismatch",float)
    myopts.add("proccol",2,"number of columns in processor grid (y direction)",int)

    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      ["booster_classic.lat",
                                       "envelope_match.cache"])
    scaling_frequency = 37.7e6
    part_per_cell = myopts.get("partpercell")
    pipe_radius = 0.04
    offset = (0,0,0)
    griddim = myopts.get("scgrid")
    proccol = myopts.get("proccol")
    if MPI.size == 1:
        proccol = 1
    num_particles = int(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell)

    ee = synergia.Error_eater()
    ee.start()
    cell_line = range(0,25)
    order = myopts.get("maporder")
    kicks = myopts.get("kickspercell")
    for cell in range(1,25):
        cell_line[cell] = Line("booster_classic.lat","bcel%02d"%cell ,0.4,
                               scaling_frequency,order,kicks,myopts)
    injcell_line = Line("booster_classic.lat","bcelinj" ,0.4,
                        scaling_frequency,order,kicks,myopts)
    part_kicks = kicks/2
    if (kicks%2 != 0):
        raise RuntimeError, "number of kicks per cell must be even"
    inja_line = Line("booster_classic.lat","bcel01a" ,0.4,
                     scaling_frequency,order,part_kicks,myopts)
    injb_line = Line("booster_classic.lat","bcel01b" ,0.4,
                     scaling_frequency,order,part_kicks,myopts)

    beam_parameters = get_beam_parameters(injcell_line,myopts)

    bunch = s2_fish.Macro_bunch(synergia.PH_NORM_mp,1)
    bunch.init_gaussian(num_particles,myopts.get("current"),beam_parameters)

    linear_map = injcell_line.gourmet.get_single_linear_map()
    correct_for_dispersion(bunch,linear_map)

    last_inj_turn = myopts.get("injturns")
    if last_inj_turn < 1:
        raise RuntimeError,"Number of injection turns must be >=1"
        
    last_turn = myopts.get("turns")

    s = 0.0
    mean = Numeric.zeros(6,'d')
    std = None
    diag = synergia.Diagnostics(injcell_line.gourmet.get_initial_u())
    diag.add(s,bunch)
    
    if myopts.get("track"):
        mytracker = synergia.Tracker('/tmp',myopts.get("trackfraction"))

    for turn in range(1,last_turn+1):
        update_rf(cell_line,myopts,turn)
        if MPI.rank==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        if turn % myopts.get("saveperiod") == 0:
            bunch.write_particles("turn_%02d.h5"%(turn-1))
        if turn == 1:
            pass
        elif (turn <= last_inj_turn):
            s = synergia.propagate(s,inja_line.gourmet,bunch,diag,griddim)
            if MPI.rank == 0:
                print "inj_a",
            sys.stdout.flush()
            inject_bunch = s2_fish.Macro_bunch(synergia.PH_NORM_mp,1)
            inject_bunch.init_gaussian(num_particles,myopts.get("current"),beam_parameters)
            correct_for_dispersion(inject_bunch,linear_map)
            bunch.inject(inject_bunch)
        else:
            s = synergia.propagate(s,cell_line[1].gourmet,bunch,diag,griddim)
            if MPI.rank == 0:
                print "01",
            sys.stdout.flush()
        if (turn<=last_inj_turn):
            s = synergia.propagate(s,injb_line.gourmet,bunch,diag,griddim)
            if MPI.rank == 0:
                print "inj_b",
            sys.stdout.flush()
        for cell in range(2,25):
            s = synergia.propagate(s,cell_line[cell].gourmet,bunch,diag,griddim)
            if MPI.rank == 0:
                print "%02d" % cell,
            sys.stdout.flush()
            if cell % 12 == 2 and myopts.get("track"):
                mytracker.add(bunch,s)
        if MPI.rank==0:
            print
    if turn % myopts.get("saveperiod") == 0:
        bunch.write_particles("turn_%02d.g5"%turn)
    diag.write_hdf5("booster_output");
    MPI.WORLD.Barrier()
    if MPI.rank == 0:
        print "elapsed time =",time.time() - t0

