#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import physics_constants
import bunch
import diagnostics
try:
    import pylab
except:
    print "pylab not available"
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

import UberPkgpy #SpaceChargePkgpy

import apply_map
import error_eater
import options

import thread

from mpi4py import MPI

from math import pi

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

        pylab.hold(0)
        pylab.ioff()

#         pylab.subplot(2,1,1)
#         pylab.plot(z,zprime,'r,')
#         pylab.title("phase space")
#         limits = pylab.axis()
#         limits[0] = -pi
#         limits[1] = pi
#         pylab.axis(limits)

        d = diagnostics.Diagnostics()

        pylab.subplot(2,2,3)
        pylab.plot(d.s/474.2,d.num_part,'g')
        pylab.title("number of particles")

        pylab.subplot(2,2,4)
        pylab.plot(d.s[1:len(d.std):4]/474.2,d.std[1:len(d.std):4,diagnostics.x])
        pylab.title("horizontal width")

        pylab.ion()
        pylab.draw()
        global plot_index
        pylab.savefig("plot%03d.png" % plot_index)
        plot_index += 1
    
rfcells = [14,15,16,17,19,21,22,23,24]

class Line:
    def __init__(self, mad_file, line_name, kinetic_energy, scaling_frequency,
                 map_order, kicks):
        global rfcells
        self.rfnames = ["bcel%02d" % cell for cell in rfcells]
        self.g = gourmet.Gourmet(mad_file,line_name, kinetic_energy,
                                 scaling_frequency, map_order)
        if line_name in self.rfnames:
            insert_rf(self.g)
        self.g.insert_space_charge_markers(2*kicks)
        if line_name in self.rfnames:
            self._get_rf_mapping()
        self.scaling_frequency = scaling_frequency
        self.kicks = kicks
        self.line_length = self.g.orbit_length()
        self.tau = 0.5*self.line_length/self.kicks
        self.name = line_name
        
    def propagate(self, s, bunch, sc_params):
        if MPI.rank == 0:
            print "%s" % self.name,
        sys.stdout.flush()
        for kick in range(0,self.kicks):
            self.g.get_fast_mapping(kick*2).apply(bunch.particles(),
                                                  bunch.num_particles_local())
            UberPkgpy.Apply_SpaceCharge_external(\
                bunch.get_beambunch(),
                sc_params.pgrid.get_pgrid2d(),
                sc_params.field.get_fieldquant(),
                sc_params.field.get_compdom(),
                sc_params.field.get_period_length(),
                sc_params.cgrid.get_bc_num(),
                sc_params.field.get_pipe_radius(),
                self.tau, 0, self.scaling_frequency)
            self.g.get_fast_mapping(kick*2+1).apply(bunch.particles(),
                                                    bunch.num_particles_local())
            s += 2.0*self.tau
            b.write_fort(s)
        return s

    def _get_rf_mapping(self):
        self.rf_strength_mapping = {}
        self.rf_phase_mapping = {}
        self.g.iterator.reset()
        element = self.g.iterator.next()
        while element:
            if element.Type() == 'thinrfcavity':
                self.rf_strength_mapping[element.Name()] = element.Strength()
                if element.Name()[0:3] == 'RFA':
                    self.rf_phase_mapping[element.Name()] = 0
                elif element.Name()[0:3] == 'RFB':
                    self.rf_phase_mapping[element.Name()] = 1
                else:
                    print "****** Warning: unknown rfcavity",element.Name()
            element = self.g.iterator.next()

    def set_rf_params(self,voltage_ev,phase_a,phase_b):
        if self.name in self.rfnames:
            self.g.iterator.reset()
            element = self.g.iterator.next()
            phases = [phase_a,phase_b]
            while element:
                if element.Type() == 'thinrfcavity':
                    strength = self.rf_strength_mapping[element.Name()] *\
                               voltage_ev
                    phase = phases[self.rf_phase_mapping[element.Name()]]
                    element.setStrength(strength)
                    element.setPhi(phase)
#                     if MPI.rank == 0:
#                         print element.Name(),
#                         print "strength =",element.Strength(),
#                         print "phi =",element.getPhi(),
#                         print "freq =",element.getFrequency()

                element = self.g.iterator.next()
            self.g.generate_mappings()
        

def get_beam_parameters(line, opts):
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

class Sc_params:
    def __init__(self, griddim, bp):
        self.pgrid = processor_grid.Processor_grid(1)
        self.cgrid = computational_grid.Computational_grid(griddim[0],
                                                      griddim[1],
                                                      griddim[2],
                                                      "trans finite, long periodic round")
        piperad = 0.04
        self.field = field.Field(bp, self.pgrid, self.cgrid, piperad)

element_keeper = []

def insert_rf(g):
    voltage = 1.0
    phase0 = 0
    g.iterator.reset()
    element = g.iterator.next()
    while element:
        if element.Name() == "LONGA:":
            rf_total_length = element.OrbitLength(g.particle)
            cavity_length = 2.56
            drift_length = (rf_total_length - 2.0*cavity_length)/2.0
            freq = 37.7e6/(2.0*math.pi)
            phase_a = phase0 - pi
            phase_b = phase0 + pi
            quality = 0
            impedance = 0
            pre_rf_drift = beamline.drift("PRERF:",drift_length)
            element_keeper.append(pre_rf_drift)
            rfda1 = beamline.drift("RFDA1:",cavity_length/2.0)
            element_keeper.append(rfda1)
            rfda2 = beamline.drift("RFDA2:",cavity_length/2.0)
            element_keeper.append(rfda2)
            rfdb1 = beamline.drift("RFDB1:",cavity_length/2.0)
            element_keeper.append(rfdb1)
            rfdb2 = beamline.drift("RFDB2:",cavity_length/2.0)
            element_keeper.append(rfdb2)
            rfa = beamline.thinrfcavity("RFA:",freq,voltage,phase_a,
                                        quality,impedance)
            element_keeper.append(rfa)
            rfb = beamline.thinrfcavity("RFB:",freq,voltage,phase_b,
                                    quality,impedance)
            element_keeper.append(rfb)
            post_rf_drift = beamline.drift("POSTRF:",drift_length)
            element_keeper.append(post_rf_drift)
            g.beamline.putAbove(element,pre_rf_drift)
            g.beamline.putAbove(element,rfda1)
            g.beamline.putAbove(element,rfa)
            g.beamline.putAbove(element,rfda2)
            g.beamline.putAbove(element,rfdb1)
            g.beamline.putAbove(element,rfb)
            g.beamline.putAbove(element,rfdb2)
            g.beamline.putAbove(element,post_rf_drift)
            element_to_remove = element
        element = g.iterator.next()
    g.beamline.remove(element_to_remove)

def update_rf(cell_line,opts,turn):
    global rfcells
    for cell in rfcells:
        if turn < (opts.get("norfturns")+opts.get("rampturns")):
            ramp = (turn-opts.get("norfturns"))/(1.0*opts.get("rampturns"))
            offset_a = -pi/2.0*(1-ramp)
            offset_b = pi/2.0*(1-ramp)
        else:
            offset_a = 0.0
            offset_b = 0.0
###        print offset_a,offset_b
        cell_line[cell].set_rf_params(opts.get("rfvoltage"),
                                      opts.get("rfphase") + offset_a,
                                      opts.get("rfphase") + offset_b)

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
    myopts.add("rfvoltage",0.6e6/18.0,"rf voltage (MV)",float)
    myopts.add("rfphase",0.0,"rf cavity phase (rad)",float)
    myopts.add("rampturns",200,"number of turns for rf phase ramping",int)
    myopts.add("norfturns",0,"number of turns without any rf BROKEN!",int)
    myopts.add("injturns",12,"myoptsnumber of injection turns",int)
    myopts.add("scgrid",[33,33,33],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)
    
    myopts.parse_argv(sys.argv)

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
                               scaling_frequency,order,kicks)
        cell_line[cell].set_rf_params(myopts.get("rfvoltage"),
                                      myopts.get("rfphase") - pi/2.0,
                                      myopts.get("rfphase") + pi/2.0)
    injcell_line = Line("booster_classic.lat","bcelinj" ,0.4,
                        scaling_frequency,order,kicks)
    part_kicks = 2
    inja_line = Line("booster_classic.lat","bcel01a" ,0.4,
                     scaling_frequency,order,part_kicks)
    injb_line = Line("booster_classic.lat","bcel01b" ,0.4,
                     scaling_frequency,order,part_kicks)

    bp = get_beam_parameters(injcell_line,myopts)

    sc_params = Sc_params(griddim,bp)
    b = bunch.Bunch(myopts.get("current"), bp, num_particles, sc_params.pgrid)
    b.generate_particles()

    #b.write_particles("initial.dat")

    last_inj_turn = myopts.get("injturns")
    last_turn = myopts.get("turns")

    s = 0.0
    b.write_fort(s)

    if MPI.rank == 0 and myopts.get("showplot"):
        pylab.ion()
    time_file = open("timing.dat","w")
    t1 = time.time()
    for turn in range(1,last_inj_turn+1):
        if MPI.rank == 0:
            time_file.write("%g %g\n" % (s,time.time()-t1))
            t1 = time.time()
            time_file.flush()
        update_rf(cell_line,myopts,turn)
        if MPI.rank==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        s = injb_line.propagate(s, b, sc_params)
        if turn % myopts.get("plotperiod") == 0 and myopts.get("showplot"):
            plot_long(b)
        for cell in range(2,25):
            s = cell_line[cell].propagate(s,b, sc_params)
            if cell % 24 == 0 and turn % myopts.get("plotperiod") == 0 \
                   and myopts.get("showplot"):
                plot_long(b)
        if not turn == last_inj_turn:
            s = inja_line.propagate(s,b,sc_params)
            binj = bunch.Bunch(myopts.get("current"), bp, num_particles,
                               sc_params.pgrid)
            binj.generate_particles()
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
        update_rf(cell_line,myopts,turn)
        if MPI.rank==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        for cell in range(1,25):
            s = cell_line[cell].propagate(s,b, sc_params)
            if cell % 24 == 0 and turn % myopts.get("plotperiod") == 0 \
                   and myopts.get("showplot"):
                plot_long(b)
        if MPI.rank==0:
            print        
        if turn % myopts.get("saveperiod") == 0:
            b.write_particles("turn_%04d.dat" % turn)
    print "elapsed time =",time.time() - t0
    time_file.close()

