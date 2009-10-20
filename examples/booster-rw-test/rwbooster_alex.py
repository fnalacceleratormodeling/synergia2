#!/usr/bin/env bwpython

import synergia
import s2_fish

import basic_toolkit
import beamline

import numpy
import numpy.linalg
import numpy
import numpy
import time
import math
import sys
from math import pi,sqrt

from mpi4py import MPI
#from fakempi import MPI
try:
    import pylab
except:
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "pylab not available"

def write_map(map):
    f = open("map.dat","w")
    for i in range(0,6):
        for j in range(0,6):
            f.write("%0.15g\n" % map[i,j]);
    f.close()
    
def read_map():
    map = numpy.zeros([6,6],'d')
    f = open("map.dat","r")
    for i in range(0,6):
        for j in range(0,6):
            map[i,j] = float(f.readline())
    f.close()
    return map
    
rfcells = [14,15,16,17,19,21,22,23,24]
#~ rfcells = []
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
                                 scaling_frequency, map_order, delay_complete=True)
        if line_name in self.rfnames:
            self.insert_rf()
        self.gourmet.insert_space_charge_markers(kicks)
        #self.gourmet.complete_setup()
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
    
def update_rf(cell_line,opts):
    global rfcells
    for cell in rfcells:
        phase0 = opts.get("rfphase")
        cell_line[cell].set_rf_params(opts.get("rfvoltage"),phase0,phase0)
      

def ha_match(map,beam_parameters,emitx,emity,dpop):
    numpy_map = map
    evals,evect_matrix = numpy.linalg.eig(numpy_map)
    evects = []
    for i in range(0,6):
        evects.append(evect_matrix[:,i])
    E = range(0,3)
    remaining = range(5,-1,-1)
    for i in range(0,3):
        # find complex conjugate among remaining eigenvectors
        first = remaining.pop()
        best = 1.0e30
        conj = -1
        for item in remaining:
            sum = evects[first]+evects[item]
            if abs(numpy.max(sum.imag)) < best:
                best = abs(numpy.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError,"failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        #~ print first,conj,best
        tmp=numpy.outer(evects[first],
            numpy.conjugate(evects[first]))
        tmp+=numpy.outer(evects[conj],
            numpy.conjugate(evects[conj]))
        E[i]=tmp.real

    for srt in range(1,3):
        if abs(numpy.linalg.det(E[srt][0:2,0:2])) > abs(numpy.linalg.det(E[0][0:2,0:2])):
            tmp = E[0]
            E[0] = E[srt]
            E[srt] = tmp
    if abs(numpy.linalg.det(E[2][2:4,2:4])) > abs(numpy.linalg.det(E[1][2:4,2:4])):
        tmp = E[1]
        E[1] = E[2]
        E[2] = tmp
    Cxy, Cxpyp, Cz, Czp = beam_parameters.get_conversions()
    gamma = beam_parameters.get_gamma()
    C = numpy.zeros([6,6],'d')
    C += E[0]*emitx/(Cxpyp*sqrt(abs(numpy.linalg.det(E[0][0:2,0:2]))))
    C += E[1]*emity/(Cxpyp*sqrt(abs(numpy.linalg.det(E[1][2:4,2:4]))))
    C += E[2]*dpop*dpop/(gamma*gamma*E[2][5,5])
    
    return C


if ( __name__ == '__main__'):
    t0 = time.time()

    myopts = synergia.Options("booster")
    myopts.add("current",0.035*12,"current",float)
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",1,"map order",int)
    myopts.add("emittance",3.0e-6,"(x-,y-)emittance",float)
    myopts.add("zfrac",0.1,"z width as fraction of bucket",float)
    myopts.add("dpop",5.0e-4,"(delta p)/p",float)
    myopts.add("kickspercell",4,"kicks per cell",int)
    myopts.add("turns",100,"number of booster revolutions",int)
    myopts.add("rfvoltage",0.6e6/18.0,"rf voltage (MV)",float)
    myopts.add("rfphase",0.0,"rf cavity phase (rad)",float)
    myopts.add("scgrid",[32,32,32],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[1,100],"fraction of particles to track (numer,denom)",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)
    myopts.add("numparticles", 0, "number of particles; if non-zero, supersedes partpercell", int)
    myopts.add("offsetx",0.00001,"x beam offset (m)",float)
    myopts.add("offsety",0.00001,"y beam offset (m)",float)
    myopts.add("latticefile","booster_classic.lat","lattice file",str)
    myopts.add("space_charge",0,"flag for space_charge",int)
    myopts.add("impedance",0,"flag for impedance",int)
   # myopts.add("pipe_symmetry","circular","",str)
    myopts.add("pipe_symmetry","x_parallel_plates","",str)
    myopts.add("pipe_radius", 0.025,"average magnet self-aperture (m)",float)
    myopts.add("pipe_conduct", 1.4e6,"conductivity # [/s] (stainless steel)",float)
    myopts.add("bunchnp",6.0e10,"number of particles per bunch",float)
    myopts.add("bunches",1,"",int)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])
    scaling_frequency = 37.77e6
    
    kinetic_energy=0.4
    offset = (0,0,0)
    griddim = myopts.get("scgrid")
    num_particles = myopts.get("numparticles")
    if num_particles == 0:
        num_particles = int(griddim[0]*griddim[1]*griddim[2] *\
                                         myopts.get("partpercell"))

    ee = synergia.Error_eater()
    ee.start()
    cell_line = range(0,25)   
    order = myopts.get("maporder")
    kicks = myopts.get("kickspercell")
    for cell in range(1,25):
        cell_line[cell] = Line(myopts.get("latticefile"),"bcel%02d"%cell ,kinetic_energy,
                               scaling_frequency,order,kicks,myopts)

        
 
    #beam_parameters = synergia.Beam_parameters(synergia.PH_NORM_mp,
                                         #1.0, 0.4, 0.0, scaling_frequency,
                                         #0)
                                         
    beam_parameters = synergia.Beam_parameters(synergia.PH_NORM_mp, 1.0, kinetic_energy,
                                         initial_phase_rad=0., scaling_frequency_Hz=scaling_frequency,
                                         transverse=0, adjust_zlength_to_freq=1)

   

    full_map = numpy.identity(7, 'd')
    update_rf(cell_line,myopts)
    for cell in range(1,25):
        cell_line[cell].gourmet.complete_setup()
    
#Line completed...    
    
    for cell in range(1,25):
        linear_map = cell_line[cell].gourmet.get_single_linear_map()
        full_map = numpy.dot(linear_map,full_map)
    
        #print "full_map ="
        #print numpy.array2string(full_map,precision=2,suppress_small=True)

    C = ha_match(full_map[0:6,0:6],beam_parameters,myopts.get("emittance"),
                 myopts.get("emittance"),myopts.get("dpop"))
    #print "C = "
    
    #print numpy.array2string(C,precision=2)
    
    # jfa: this is an ugly hack
    beam_parameters.offset_x_m = myopts.get("offsetx")
    beam_parameters.offset_y_m = myopts.get("offsety")
#    bunch.init_gaussian_covariance(num_particles, myopts.get("current"), beam_parameters, C)
    
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "gamma=",beam_parameters.get_gamma(),"  beta=", beam_parameters.get_beta(),"  betagamma=",betagamma
    
    #bunchnp=myopts.get("current")/(synergia.physics_constants.PH_MKS_e *beam_parameters.scaling_frequency_Hz)
    #print "bunch_np=",bunchnp
    
    
    numbunches = myopts.get("bunches")
    bunchnp=myopts.get("bunchnp")
   # covariance=numpy.zeros([6,6],'d')
    bunches = []
    
    for bunchnum in range(0,numbunches):
        diag=synergia.Diagnostics(cell_line[1].gourmet.get_initial_u())
        covariance=C.copy()
        bunches.append(s2_fish.Macro_bunch.gaussian_covariance(bunchnp,num_particles,beam_parameters,covariance,diagnostics=diag,bucket_num=bunchnum,periodic=True))        
        bunches[bunchnum].write_particles("begin-%02d"%bunchnum)
#************** comment the following*************   
    #pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    #beam_parameters.x_params(sigma = 0.0002, lam = 0.008 * pz,
                             #r = 0.5,offset=myopts.get("offsetx"), offset_p =0.)
    #beam_parameters.y_params(sigma = 0.001, lam = 0.001* pz,
                             #r = 0.5,offset=myopts.get("offsety"), offset_p = 0. )
    #beam_parameters.z_params(sigma = 0.7, 
                             #lam = 0.0002* pz,  offset=0.,
                             #offset_p = 0.)                                                 
    #bunch= s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,periodic=True)
    
 #**************************************************** 
   # bunch_sp=2.0*math.pi*beam_parameters*synergia.physics_constants.PH_MKS_c/beam_parameters.get_omega()
    bunch_sp=beam_parameters.get_z_length()  
    mbunches=s2_fish.Multiple_bunches(bunches, bunch_sp) 
   
    
    pipe_radius = myopts.get("pipe_radius")
    space_charge=myopts.get("space_charge")
    if space_charge:
        griddim = myopts.get("scgrid")
        solver="s2_fish_cylindrical"
       # griddim = (16,16,33)
       # solver="s2_fish_3d"

        sp_ch=s2_fish.SpaceCharge(solver,griddim,radius_cylindrical=pipe_radius,periodic=True) 
        if MPI.COMM_WORLD.Get_rank() ==0: 
            print " sp_ch grid=",sp_ch.get_grid()
            print " sp_ch solver=",sp_ch.get_solver()
            print " sp_ch pipe radius=",sp_ch.get_radius_cylindrical()
    else:
       sp_ch=None       
    
    impedance=myopts.get("impedance")
    if impedance:
        pipe_conduct= myopts.get("pipe_conduct") # [ohm^-1 m^-1] (stainless steel)
        prev_turns=10  
        wall_thickness=0.114        
        pipe_symmetry=myopts.get("pipe_symmetry")
        kick="full"
        lgridnum=40
        line_length =0.
        
        for cell in range(1,25):
            line_length += cell_line[cell].gourmet.orbit_length() 
        
        rw_impedance=s2_fish.Impedance(pipe_radius, pipe_conduct,wall_thickness, line_length,lgridnum,
             pipe_symmetry=pipe_symmetry,paking_frac=0.6,kick=kick,nstored_turns=prev_turns)
                #pipe_symmetry="x_parallel_plates") 
        if MPI.COMM_WORLD.Get_rank() ==0: 
            print "IMPEDANCE PIPE radius=", rw_impedance.get_pipe_radius()
            print "IMPEDANCE PIPE wall_thickness=",rw_impedance.get_wall_thickness()
            print "IMPEDANCE PIPE symmetry=",rw_impedance.get_pipe_symmetry()
            print " Orbith length=",rw_impedance.get_orbit_length()
    else:
        rw_impedance=None     
    
    
    log = open("log","w")
    if MPI.COMM_WORLD.Get_rank() ==0:
            output = "start propagation"
            print output
            log.write("%s\n" % output)
            log.flush() 
   
    s = 0.0
   # bunch.add_diagnostics(s)

  
        
    if myopts.get("track"):
        mytracker = synergia.Tracker('/tmp',myopts.get("trackfraction"))
        mytracker.add(bunch,s)

    for turn in range(1,myopts.get("turns")+1):
        t1 = time.time()
        #~ dt0 = time.time()
        #~ diag.write_hdf5("booster_output")
        #~ dt1 = time.time()
        #~ print "writing booster_output took",dt1-dt0,"s"
        if MPI.COMM_WORLD.Get_rank()==0:
            print "turn %d:" % turn,
            sys.stdout.flush()
        if turn % myopts.get("saveperiod") == 0:
            for bunchnum in range(0,numbunches):
                bunches[bunchnum].write_particles("bunch%02d_turn_%03d.h5" %(bunchnum, (turn-1)))
        for cell in range(1,25):
            
            s=synergia.propagate(s,cell_line[cell].gourmet, mbunches, space_charge=sp_ch,impedance=rw_impedance)
          
            if MPI.COMM_WORLD.Get_rank() == 0:
                print "%02d" % cell,
            sys.stdout.flush()
            if cell % 12 == 0 and myopts.get("track"):
                mytracker.add(bunch,s)
        if MPI.COMM_WORLD.Get_rank() ==0:
            print
            output = "turn %d time = %g"%(turn,time.time() - t1)    
            print output
            log.write("%s\n" % output)
            log.flush()   
            print
            
    #if turn % myopts.get("saveperiod") == 0:
        #bunch.write_particles("turn_%03d.g5"%turn)
    for bunchnum in range(0,numbunches):
        bunches[bunchnum].diagnostics.write_hdf5("booster_output-%02d"%bunchnum)
    for bunchnum in range(0,numbunches):
        bunches[bunchnum].write_particles("end-%02d"%bunchnum)   
    
    
    #if myopts.get("track"):
        #mytracker.close()
        #mytracker.show_statistics() 
   # print "      BEFORE barrier on rank= ", MPI.COMM_WORLD.Get_rank()     
    MPI.WORLD.Barrier()
   # print " after barrier on rank= ", MPI.COMM_WORLD.Get_rank()
    log.close()
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "elapsed time =",time.time() - t0
    
