#!/usr/bin/env bwpython

import synergia
import s2_fish

import basic_toolkit
import beamline

import Numeric
import LinearAlgebra
import numpy
import MLab
import time
import math
import sys
from math import pi,sqrt

from s2_fish.s2_containers import Real_scalar_field
from s2_fish.s2_deposit import rho_to_rwvars
from s2_fish.s2_deposit import deposit_charge_cic

#from mpi4py import MPI
from fakempi import MPI
try:
    import pylab
except:
    if MPI.rank == 0:
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
                best = abs(MLab.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError,"failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        #~ print first,conj,best
        tmp=Numeric.outerproduct(evects[first],
            Numeric.conjugate(evects[first]))
        tmp+=Numeric.outerproduct(evects[conj],
            Numeric.conjugate(evects[conj]))
        E[i]=tmp.real

    for srt in range(1,3):
        if abs(LinearAlgebra.determinant(E[srt][0:2,0:2])) > abs(LinearAlgebra.determinant(E[0][0:2,0:2])):
            tmp = E[0]
            E[0] = E[srt]
            E[srt] = tmp
    if abs(LinearAlgebra.determinant(E[2][2:4,2:4])) > abs(LinearAlgebra.determinant(E[1][2:4,2:4])):
        tmp = E[1]
        E[1] = E[2]
        E[2] = tmp
    Cxy, Cxpyp, Cz, Czp = beam_parameters.get_conversions()
    gamma = beam_parameters.get_gamma()
    C = Numeric.zeros([6,6],'d')
    C += E[0]*emitx/(Cxpyp*sqrt(abs(LinearAlgebra.determinant(E[0][0:2,0:2]))))
    C += E[1]*emity/(Cxpyp*sqrt(abs(LinearAlgebra.determinant(E[1][2:4,2:4]))))
    C += E[2]*dpop*dpop/(gamma*gamma*E[2][5,5])
    
    return C


if ( __name__ == '__main__'):
    t0 = time.time()

    myopts = synergia.Options("booster")
    myopts.add("current",0.035*12,"current",float)
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",3.0e-6,"(x-,y-)emittance",float)
    myopts.add("zfrac",0.1,"z width as fraction of bucket",float)
    myopts.add("dpop",5.0e-4,"(delta p)/p",float)
    myopts.add("kickspercell",4,"kicks per cell",int)
    myopts.add("turns",100,"number of booster revolutions",int)
    myopts.add("rfvoltage",0.6e6/18.0,"rf voltage (MV)",float)
    myopts.add("rfphase",0.0,"rf cavity phase (rad)",float)
    myopts.add("scgrid",[33,33,33],"Space charge grid",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("track",0,"whether to track particles",int)
    myopts.add("trackfraction",[1,100],"fraction of particles to track (numer,denom)",int)
    myopts.add("partpercell",1.0,"average number of particles per cell",float)
    myopts.add("numparticles", 0, "number of particles; if non-zero, supersedes partpercell", int)
    myopts.add("offsetx",0,"x beam offset (m)",float)
    myopts.add("offsety",0,"y beam offset (m)",float)
    myopts.add("latticefile","booster_classic.lat","lattice file",str)
    myopts.add("space_charge",0,"flag for space_charge",int)
    myopts.add("impedance",1,"flag for impedance",int)
    myopts.add("pipe_conduct", 1.4e6,"conductivity # [/s] (stainless steel)",float)

    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])
    scaling_frequency = 37.7e6
    pipe_radius = 0.04
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
        cell_line[cell] = Line(myopts.get("latticefile"),"bcel%02d"%cell ,0.4,
                               scaling_frequency,order,kicks,myopts)

    beam_parameters = synergia.Beam_parameters(synergia.PH_NORM_mp,
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         0)

    bunch = s2_fish.Macro_bunch(synergia.PH_NORM_mp,1)
#    bunch.init_gaussian(num_particles,myopts.get("current"),beam_parameters)
    full_map = Numeric.identity(7, 'd')
    update_rf(cell_line,myopts)
    for cell in range(1,25):
        linear_map = cell_line[cell].gourmet.get_single_linear_map()
        full_map = Numeric.matrixmultiply(linear_map,full_map)
    
    #~ print "full_map ="
    #~ print Numeric.array2string(full_map,precision=2,suppress_small=True)

    C = ha_match(full_map[0:6,0:6],beam_parameters,myopts.get("emittance"),
                 myopts.get("emittance"),myopts.get("dpop"))
    #~ print "C = "
    #~ print Numeric.array2string(C,precision=2)

    # jfa: this is an ugly hack
    beam_parameters.offset_x_m = myopts.get("offsetx")
    beam_parameters.offset_y_m = myopts.get("offsety")
    bunch.init_gaussian_covariance(num_particles, myopts.get("current"), beam_parameters, C)

    s = 0.0
    diag = synergia.Diagnostics(cell_line[1].gourmet.get_initial_u())
    diag.add(s,bunch)
    
    
    print "num_particles =", bunch.get_num_particles_local()
    print "means =",diag.means[0]
    
    shape = griddim
    means, stds = synergia.get_spatial_means_stds(bunch)
    n_sigma = 8.0
    size = list(n_sigma*stds)
    offset = list(means)
    periodic = False
    print "shape =",shape
    print "size =",size
    print "offset =",offset
    rho = Real_scalar_field(shape,size,offset)
    deposit_charge_cic(rho,bunch.get_store(),periodic)
    zdensity = Numeric.zeros((shape[2]),'d')
    xmom = Numeric.zeros((shape[2]),'d')
    ymom= Numeric.zeros((shape[2]),'d')
    rho_to_rwvars(rho,zdensity,xmom,ymom)

    #~ rho.get_points().print_("rho")

    print "zdensity =",zdensity
    print "sum zdensity =",sum(zdensity)
    cell_size = rho.get_cell_size()
    print "cell_size =",cell_size
    cell_volume = cell_size[0]*cell_size[1]*cell_size[2]
    print "cell_volume =",cell_volume
    print "xmom =",xmom
    print "ymom =",ymom
