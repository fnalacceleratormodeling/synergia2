#!/usr/bin/env python
import sys
import synergia
from circular_options import opts
import numpy as np
from synergia.optics.one_turn_map import linear_one_turn_map
from mpi4py import MPI

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        print "error, map is unstable"
    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo_s.lat")
lattice_length = lattice.get_length()
reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
freq = opts.harmno * beta * synergia.foundation.pconstants.c/lattice_length

if MPI.COMM_WORLD.Get_rank() ==0:
    print "lattice length: ", lattice_length
    print "reference particle energy: ", energy
    print "reference particle beta: ", beta
    print "reference particle gamma: ", gamma
    print "RF freq: ", freq
# set rf cavity frequency
# harmno * beta * c/ring_length








# Don't need that?  Lattice has voltage set
# rf cavity voltage, 
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", opts.rf_voltage)
        elem.set_double_attribute("freq", freq)
       # print" atributes=", elem.get_double_attributes()
        #elem.set_double_attribute("lag", 0.5)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
map = linear_one_turn_map(lattice_simulator)
if MPI.COMM_WORLD.Get_rank() ==0:
    print "bucket length=",lattice_simulator.get_bucket_length()
    print "num of buckets=",lattice_simulator.get_number_buckets()
    print "Linear one turn map"
    print np.array2string(map,max_line_width=200)






l,v = np.linalg.eig(map)
if MPI.COMM_WORLD.Get_rank() ==0:
    print "eigenvalues of one turn map: ", l
    print "absolute values of eigenvalues (should all be 1): ", abs(l)
    print "fractional tunes from eigenvalues: ", np.log(l).imag/(2.0*np.pi)


[[ax,ay],[bx,by]] = synergia.optics.get_alpha_beta(map)




if MPI.COMM_WORLD.Get_rank() ==0:
    print "Lattice functions assuming uncoupled map:"
    print "alpha x: ", ax
    print "alpha y: ", ay
    print "beta x: ", bx
    print "beta y: ", by

[az, bz, qz] = map2twiss(map[4:6,4:6])
if MPI.COMM_WORLD.Get_rank() ==0:
    print "alpha z (better be small): ", az
    print "beta z: ", bz





emit = opts.norm_emit
if MPI.COMM_WORLD.Get_rank() ==0:
    print "generating particles with transverse emittance: ", emit

rms_index=[0,2,4]
arms=np.sqrt(emit*bx)
brms=np.sqrt(emit*by)
crms=opts.stdz

#rms_index=[1,3,5]
#arms=0.002691405845277
#brms=0.000691411827683
#crms=4.84505679134e-06


covar = synergia.optics.matching._get_correlation_matrix(map,arms,brms,crms,beta,rms_index)
if MPI.COMM_WORLD.Get_rank() ==0:
   # print "covariance matrix"
   # print np.array2string(covar,max_line_width=200)
    print "stdx =",np.sqrt(emit*bx)," stdy= ", np.sqrt(emit*by)

                                               
#bunch= synergia.optics.generate_matched_bunch(lattice_simulator,
                                               #arms,brms,crms,
                                               #opts.num_real_particles,
                                               #opts.num_macro_particles,rms_index,
                                               #seed=opts.seed)





 #for bunchnum in range(0,numbunches):
        #diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True)
        #bunchnp=bunchnp0#*(bunchnum+1)*0.5 # bucket_num =2 in front of bucket_num =3
        #bunches.append(s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,diagnostics=diag,bucket_num=2*bunchnum,periodic=True))
        #bunches[bunchnum].write_particles("begin-%02d"%bunchnum)
        #print " bunch(",bunchnum,") periodicity=",bunches[bunchnum].periodic
       ## print "  initial means bunch(",bunchnum,")=",numpy.array(bunches[bunchnum].diagnostics.get_means())

    #print " **********************************************************************"  
    #mbunches=s2_fish.Multiple_bunches(bunches, bunch_sp)
        

bunchsp=lattice_simulator.get_bucket_length()
num_bunches=emit = opts.num_bunches



#bunch_diag_train=synergia.bunch.Bunch_with_diagnostics_train(num_bunches,bunchsp, MPI.COMM_WORLD)
#for bunchnum in range(0,num_bunches):
    #if bunch_diag_train.is_on_this_rank(bunchnum):
        #commx=bunch_diag_train.get_comm(bunchnum)
        #bunch= synergia.optics.generate_matched_bunch(lattice_simulator,
                                                #arms,brms,crms,
                                                #opts.num_real_particles*(bunchnum+1),
                                                #opts.num_macro_particles,rms_index,
                                                #seed=opts.seed, bunch_index=bunchnum,comm=commx, periodic=True)
        #particles = bunch.get_local_particles()
        ## apply offset to bunch
        #particles[:,0] = particles[:,0]+opts.x_offset
        #particles[:,2] = particles[:,2]+opts.y_offset
        #particles[:,4] = particles[:,4]+opts.z_offset     
        #diagnostics_writer_step = synergia.bunch.Diagnostics_full2(bunch, "circular_full2-%02d.h5"%bunchnum) 
        #diagnostics_writer_turn = synergia.bunch.Diagnostics_particles(bunch,"circular_particles-%02d.h5"%bunchnum,0,0,100)
        #bunch_diag=synergia.bunch.Bunch_with_diagnostics(bunch,diagnostics_writer_step, diagnostics_writer_turn)
        #bunch_diag_train.set_bunch_diag_sptr(bunchnum, bunch_diag)
        #real_num=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_real_num()
        #bucket_index=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_bucket_index()
        #if commx.Get_rank() ==0:
            #print "bunch # ",bunchnum ,"  number of particles= ",real_num, " bucket =", bucket_index
#if MPI.COMM_WORLD.Get_rank() ==0:
    #print "train bunch space=",bunch_diag_train.get_bunch_separation()


bunch= synergia.optics.generate_matched_bunch(lattice_simulator,
                                                arms,brms,crms,
                                                opts.num_real_particles,
                                                opts.num_macro_particles,rms_index,
                                                seed=opts.seed, periodic=True)
diagnostics_actions = synergia.simulation.Standard_diagnostics_actions()
#diagnostics_actions.add_per_step(synergia.bunch.Diagnostics_full2(bunch, "step_full2a.h5"))
#diagnostics_actions.add_per_turn(synergia.bunch.Diagnostics_particles(bunch, "turn_particles.h5a",0,0,100))
#bunch_with_diag=synergia.bunch.Bunch_with_diagnostics(bunch, diagnostics_actions)
#bunch_with_diag.check_bunch_pointer_in_diagnostics() 


bunch_with_diag=synergia.bunch.Bunch_with_diagnostics(bunch, diagnostics_actions)
bunch_with_diag.add_per_step_diagnostics(synergia.bunch.Diagnostics_full2(bunch, "step_full2b.h5"))
bunch_with_diag.add_per_turn_diagnostics(synergia.bunch.Diagnostics_particles(bunch, "turn_particlesb.h5",0,0,100))

diagnostics_writer_step = synergia.bunch.Diagnostics_full2(bunch, "circular_full2.h5") 
diagnostics_writer_turn = synergia.bunch.Diagnostics_particles(bunch,"circular_particles.h5",0,0,100)




no_op = synergia.simulation.Dummy_collective_operator("stub")
zgrid=40
imped= synergia.collective.Impedance("BoosterF_wake.dat",lattice_length, bunchsp,zgrid, "circular",60)
#imped= synergia.simulation.Dummy_collective_operator("stub")
impedance=opts.impedance
if impedance:
    stepper = synergia.simulation.Split_operator_stepper(
                            lattice_simulator, imped, opts.num_steps)
else:
    stepper = synergia.simulation. Split_operator_stepper(
                            lattice_simulator, no_op, opts.num_steps)                           
    #stepper = synergia.simulation.Independent_stepper_elements(


if MPI.COMM_WORLD.Get_rank() ==4:
    print "expect std_x: ", np.sqrt(emit*bx)
    print "generated std_x: ", np.std(particles[:,0])
    print "expect std_y: ", np.sqrt(emit*by)
    print "generated std_y: ", np.std(particles[:,2]);
    print "expected std_z: ", opts.stdz
    print "generated std_z: ", np.std(particles[:,4])
    print "expected std(dpop): ", opts.stdz/bz
    print "generated std(dpop): ", np.std(particles[:,5])

#diagnostics_writer_step = synergia.bunch.Diagnostics_full2(bunch, "circular_full2.h5")

#diagnostics_writer_turn = synergia.bunch.Diagnostics_particles(bunch,"circular_particles.h5",0,0,100)



#
#if MPI.COMM_WORLD.Get_rank() ==0:
#    print " real num in bunch&diag=",bunch_with_diag.get_bunch_sptr().get_real_num()




propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch_with_diag, opts.num_turns, opts.verbose)
#propagator.propagate(bunch_diag_train, opts.num_turns, opts.verbose)
#propagator.propagate(bunch, opts.num_turns, diagnostics_writer_step, diagnostics_writer_turn, opts.verbose)
#propagator.propagate(bunch, opts.num_turns, diagnostics_actions, opts.verbose)
