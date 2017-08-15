#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import synergia
from iota_options import opts
import numpy as np
import string
from synergia.optics.one_turn_map import linear_one_turn_map
from synergia.utils import Commxx
from synergia.bunch import Core_diagnostics
from mpi4py import MPI



lattices = {}
lattice_repo = './'
lattices['t1_1IO_66'] = lattice_repo + "lattice_1IO_center.madx" #centered t1 6.6 1IO lattice
lattices['t3_1IO_66'] = lattice_repo + "lattice_1IO_nll_center.madx" #centered t3 6.6 1IO lattice


if opts.nnl:#Flag for using the version of the lattice with the nonlinear element
   lname = 't3_1IO_66'
else:
   lname = 't1_1IO_66' 

lattice = synergia.lattice.MadX_reader().get_lattice("iota", lattices[lname])  

if MPI.COMM_WORLD.Get_rank() ==0:
     print "lattice from the file: ",  lattices[lname], "  made"
     


for elem in lattice.get_elements():  
    element_name=elem.get_name()
    element_type=elem.get_type()
    
    if opts.chef_propagate:
          elem.set_string_attribute("extractor_type", "chef_propagate")
    elif  opts.chef_map: 
          elem.set_string_attribute("extractor_type", "chef_map")
          
    if opts.if_aperture:
        name1=element_name[0]
        elem.set_string_attribute("aperture_loss","aperture_loss_file");  
        if (element_type=="sbend") and (name1=="m"):
                #print "name=",element_name, " h apert=",opts.get("aperture_bending")
                elem.set_string_attribute("aperture_type","rectangular")
                elem.set_double_attribute("rectangular_aperture_width",  2*opts.get("aperture_bending"))
                elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("aperture_bending"))
        else:  
                #print "name=",element_name, " h apert=",opts.get("aperture_straight")           
                elem.set_string_attribute("aperture_type","circular")
                elem.set_double_attribute("circular_aperture_radius", opts.get("aperture_straight"))


reference_particle = lattice.get_reference_particle() 
lattice_length=lattice.get_length()    
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
energy = reference_particle.get_total_energy()
if MPI.COMM_WORLD.Get_rank() ==0:
   print "  ***     before stepper   ***     "
   print " beta=",beta
   print "gamma=",gamma
   print "lattice length=",lattice_length
   print "energy=",energy
   print "  ***      ***********     ***     "  
    

#*****************************************************************************************

  

#make collective operators
#space charge
if opts.space_charge_rec:   
     grid_shape_straight=opts.scgrid_straight[:]
     radius_straight=opts.aperture_straight;
     pipe_size_straight=[2*radius_straight,2*radius_straight,lattice_length/opts.harmon]
     spcr_straight=synergia.collective.Space_charge_rectangular(pipe_size_straight, grid_shape_straight)
     if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "SPACE CHARGE SOLVER FOR RECTANGULAR PIPES"
        print "pipe_size straight section=",spcr_straight.get_pipe_size()
        print "grid for spc rectangular straight section=",spcr_straight.get_grid_shape() 
        print "________________________________________________"
        
     grid_shape_bending=opts.scgrid_bending[:]
     radius_bending=opts.aperture_bending;
     pipe_size_bending=[2*radius_bending,2*radius_bending,lattice_length/opts.harmon]
     spcr_bending=synergia.collective.Space_charge_rectangular(pipe_size_bending, grid_shape_bending)
     if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "SPACE CHARGE SOLVER FOR RECTANGULAR PIPES"
        print "pipe_size bending section=",spcr_bending.get_pipe_size()
        print "grid for spc rectangular bending section=",spcr_bending.get_grid_shape() 
        print "________________________________________________"     

elif opts.space_charge_3dh:
     grid_shape=opts.scgrid_3dh[:]
     longitudinal_kicks=True
     periodic_z = False
     z_period=gamma*lattice_length/opts.harmon
     grid_entire_period = False
     nsigma=8.
     commxx_divider= synergia.utils.Commxx_divider(opts.spc_comm_size, False)
     spc3dh=synergia.collective.Space_charge_3d_open_hockney(commxx_divider,grid_shape,longitudinal_kicks,periodic_z,z_period,grid_entire_period,nsigma )
     if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "3D HOCKNEY SOLVER WITH OPEN BOUNDARIES"      
        print "grid for spc 3dh=",grid_shape 
        print "longitudinal kicks =",longitudinal_kicks
        print "longitudinal periodicity=",periodic_z
        print "nsigma=",nsigma                    
        print "________________________________________________"
        
elif opts.space_charge_2dh: 
     grid_shape=opts.scgrid_3dh[:]
     commxx_divider= synergia.utils.Commxx_divider(opts.spc_comm_size, False)
     spc2dh=synergia.collective.Space_charge_2d_open_hockney(commxx_divider,grid_shape)
     if MPI.COMM_WORLD.Get_rank() ==0:        
        print
        print "2D HOCKNEY SOLVER WITH OPEN BOUNDARIES"      
        print "grid for spc 2dh=",grid_shape 
        print "________________________________________________"
     
     
#impedance
if opts.impedance:
    wn=opts.wave_number[:] # should be [0,0,0] unless multibunch simulations are considered
    imped_straight=synergia.collective.Impedance(opts.wakefile_straight, opts.waketype, opts.zgrid, lattice_length, float(lattice_length/opts.harmon), opts.registred_turns,opts.full_machine,wn)
    if MPI.COMM_WORLD.Get_rank() ==0:  
          print
          print "WAKES FOR STRAIGHT read from ",imped_straight.get_wake_field().get_wake_file_name()
          print "STRAIGHT orbith length=",imped_straight.get_orbit_length()
          print "STRAIGHT z_grid=",imped_straight.get_z_grid()
          print "STRAIGHT stored turns=",imped_straight.get_nstored_turns()
          print "________________________________________________"

dummy_operator=synergia.simulation.Dummy_collective_operator("dummy_impedance") # this operator  does nothing
  
     
# make the stepper
# when collective operators are present we choose  opts.steps_per_sbend kicks on each bending magnet and  opts.num_steps_straight kicks
# evenly distributed in space (s) elswhere. This is done by defining a list_choice_map.
# example: if on the element named "bobo" you desire to put 4 kicks (spaced evenly in s), and each kick has three collective operatrs, op1,op2 and op3
# you  make kicks_bobo=Kicks([op1,op2,op3],4) and list_choice_map["bobo"]=kicks_bobo.
# When the stepper is created, 4 kicks are put on the element with the name "bobo".
# After finishing with all the special elements (like bobo), you want to put 10 (evenly in s) kicks containg op4 and op5  
# over the remaining part of the ring
# you make kicks_remaining=Kicks[op4,op5],10) and list_choice_map["else"]=kicks_remaining
# "else" is a keyword for whatever reamains in the lattice, the stepper constructor search for it.
if (opts.impedance) or  (opts.space_charge_rec) or   (opts.space_charge_3dh)  or   (opts.space_charge_2dh) :
    list_choice_map=synergia.simulation.List_choice_map()
    operators_bending=[] # this contains the collective operators on the bending magnets
    if opts.impedance: operators_bending.append(dummy_operator) # no impedance on bending sector 
    if opts.space_charge_rec: operators_bending.append(spcr_bending)
    if opts.space_charge_3dh: operators_bending.append(spc3dh)
    if opts.space_charge_2dh: operators_bending.append(spc2dh)
    kicks_bending=synergia.simulation.Kicks(operators_bending,opts.steps_per_sbend)
    for elem in lattice.get_elements():  
        element_name=elem.get_name()
        element_type=elem.get_type()
        if (element_type=="sbend") and (element_name[0]=="m"):
            #print "element_name=",element_name
            list_choice_map[element_name]=kicks_bending
    
    operators_straight=[]# this contains the collective operators on the straight section
    if opts.impedance: operators_straight.append(imped_straight)
    if opts.space_charge_rec: operators_straight.append(spcr_straight)
    if opts.space_charge_3dh: operators_straight.append(spc3dh)
    if opts.space_charge_2dh: operators_straight.append(spc2dh)                                                    
    kicks_straight=synergia.simulation.Kicks(operators_straight, opts.num_steps_straight)
    list_choice_map["else"]=kicks_straight
    stepper=synergia.simulation.Split_operator_stepper_choice(lattice, opts.map_order, list_choice_map, True)
    if MPI.COMM_WORLD.Get_rank() ==0: 
       print
       print "Split_operator_stepper_choice created"             
       print "steps_per_bending magnet=",opts.steps_per_sbend
       print  "num_per_straight section=",opts.num_steps_straight   
       print  "___________________________________________________________"
else:
    stepper=synergia.simulation.Independent_stepper(lattice, opts.map_order, opts.num_steps)
    if MPI.COMM_WORLD.Get_rank() ==0: 
       print
       print  "no collective effects, no bpm, independent stepper"
       print  "___________________________________________________________" 
    
 
 
stepper.get_lattice_simulator().register_closed_orbit() 
stepper.get_lattice_simulator().set_rf_bucket_length()
stepper.get_lattice_simulator().print_lattice_functions()  # make sure if have rf volatge zero for meaningful lattice functions (in case you have rf)
if MPI.COMM_WORLD.Get_rank() ==0:   
     print "stepper:lattice_simulator: map order=",stepper.get_lattice_simulator().get_map_order()
     print  "stepper: number of steps=",len(stepper.get_steps())
     print  "stepper rf frequency=",stepper.get_lattice_simulator().get_rf_frequency()
    # stepper.print_();


reference_particle=stepper.get_lattice_simulator().get_lattice().get_reference_particle()
lattice_length=stepper.get_lattice_simulator().get_lattice().get_length()       
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
energy = reference_particle.get_total_energy()
clo=stepper.get_lattice_simulator().get_closed_orbit()

if MPI.COMM_WORLD.Get_rank() ==0: 
   print "    ********    stepper  lattice  ************     "
   print
   print " closed orbit= ",clo
   print " beta=", beta 
   print " gamma=", gamma 
   print " lattice_length=", lattice_length 
   print " closed_orbit_length=", stepper.get_lattice_simulator().get_closed_orbit_length() 
   print " energy=", energy, "GeV"   
   print " reference momentum=", reference_particle.get_momentum(), " GeV/c"     
   print 
   print "    ***********************************     "
   print
     
if opts.map_order==1:
   one_turn_map=stepper.get_lattice_simulator().get_linear_one_turn_map()
   eigenvlaues_map, eigenvectors_map=np.linalg.eig(one_turn_map)
   # [[0,2,4] when calling the get_correlation_matrix coresponds to xrms,yrms ans zrms. [1, 3, 5] would be xprms, yprms, zprms!
   correlation_matrix=synergia.bunch.get_correlation_matrix(one_turn_map,opts.xrms, opts.yrms, opts.zrms, beta, [0,2,4]); 
   if MPI.COMM_WORLD.Get_rank() ==0: 
      # print "correlation matrix for the input rms values:=",correlation_matrix
       print
       print "_______________________________________________________"
       print "one turn map=",one_turn_map
       print 
       print "absolute eigenvlaues of one turn map=",np.absolute(eigenvlaues_map) # should all be one
       print "fractional tunes=",np.absolute(np.imag(np.log(eigenvlaues_map))/(2*np.pi))
       print "1-fractional tunes=",1.-np.absolute(np.imag(np.log(eigenvlaues_map))/(2*np.pi))
       print "_______________________________________________________"
if opts.tunes_and_chroms:
    chef_frac_tunex=stepper.get_lattice_simulator().get_horizontal_tune()
    chef_frac_tuney=stepper.get_lattice_simulator().get_vertical_tune()
    chef_eigen_tunex=stepper.get_lattice_simulator().get_horizontal_tune(True)
    chef_eigen_tuney=stepper.get_lattice_simulator().get_vertical_tune(True)
    horizontal_chromaticity=stepper.get_lattice_simulator().get_horizontal_chromaticity()
    vertical_chromaticity=stepper.get_lattice_simulator().get_vertical_chromaticity()
    momentum_compaction=stepper.get_lattice_simulator().get_momentum_compaction()
    slip_factor=stepper.get_lattice_simulator().get_slip_factor()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "chef FracTune x: ", chef_frac_tunex, ", EigenTune x: ", chef_eigen_tunex
        print "chef FracTune y: ", chef_frac_tuney, ", EigenTune y: ", chef_eigen_tuney
        print " horizontal chromaticity: ", horizontal_chromaticity
        print " vertical   chromaticity: ", vertical_chromaticity
        print " momentum compaction: ", momentum_compaction
        print " slip factor: ", slip_factor  
        print 
    
       
#make bunches 
num_bunches=opts.num_bunches
bunches = []
parent_comm = Commxx()
comms = synergia.utils.generate_subcomms(parent_comm, opts.num_bunches) 
for i in range(0, num_bunches):
        commx=comms[i]
        bunch=synergia.bunch.Bunch(stepper.get_lattice_simulator().get_lattice().get_reference_particle(),
                                       opts.num_macroparticles, opts.num_real_particles, commx);
        bunch.set_bucket_index(i)
        bunch.set_z_period_length(lattice.get_length()/opts.harmon)   
        if commx.has_this_rank():
                if opts.load_bunch:              
                    bunch.read_file(opts.input_bunch_h5file);                
                    if commx.get_rank()==0: 
                        print " bunch ",i," loaded from:  ","iota_input_particles.h5"
                else:            
                    dist=synergia.foundation.Random_distribution(opts.seed,commx)
                    input_means=np.zeros(6)
                    if opts.map_order==1:
                        synergia.bunch.populate_6d(dist, bunch, input_means, correlation_matrix)
                    else:
                        raise RuntimeError,  "iota.py: bunch populate for high order maps nor implemented yet " 

                if opts.coasting_beam:   
                    dist=synergia.foundation.Random_distribution(opts.seed+13,commx)
                    synergia.bunch.populate_longitudinal_uniform(dist, bunch, bunch.get_z_period_length()) # z_period_length shoule be lattice_length
                                        
        bunches.append(bunch)       
                
                
#make bunch train      
bunch_train=synergia.bunch.Bunch_train(bunches, stepper.get_lattice_simulator().get_bucket_length());                
bunch_train_simulator = synergia.simulation.Bunch_train_simulator(bunch_train)                

# put diagnostics for bunches, make offsets and print emittances and othwer beam properties 
if opts.spc_tuneshift: 
       stepper.cs_step_lattice_functions()# required to calculate beta functions for each step
       stepper.print_cs_step_betas() 

for i in range(0, bunch_train.get_size()):       
      bunch_train_simulator.add_per_step(i, synergia.bunch.Diagnostics_full2("step_full2_%d.h5" % i))  
      bunch_train_simulator.add_per_turn(i, synergia.bunch.Diagnostics_particles("turn_particles_%d.h5" % i), opts.turn_period)
      if opts.bulk_track:
           bunch_train_simulator.add_per_turn(i,synergia.bunch.Diagnostics_bulk_track("bulk_track_%d.h5" % i, opts.num_macroparticles))
      if opts.space_charge_rec:
           if bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm().has_this_rank(): 
                comm_spc= synergia.utils.make_optimal_spc_comm(bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm(),opts.spc_comm_size,0)
                spcr_bending.set_fftw_helper(comm_spc,0)
                spcr_straight.set_fftw_helper(comm_spc,0)
                if opts.spc_tuneshift: 
                     sp_diagnostics=synergia.collective.Diagnostics_space_charge_rectangular("space_charge_rec_diagnostics_%d.h5" % i)
                     sp_diagnostics.set_bunch(bunch_train_simulator.get_bunch_train().get_bunches()[i])                  
                     spcr_bending.add_diagnostics(sp_diagnostics)
                     spcr_straight.add_diagnostics(sp_diagnostics)
      
      elif opts.space_charge_3dh:
          if opts.spc_tuneshift: 
               sp_diagnostics=synergia.collective.Diagnostics_space_charge_3d_hockney("space_charge_3dh_diagnostics_%d.h5" % i)
               sp_diagnostics.set_bunch(bunch_train_simulator.get_bunch_train().get_bunches()[i]) 
               spc3dh.add_diagnostics(sp_diagnostics);  
     
      if opts.apertures_loss:
          if bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm().has_this_rank(): 
               diag_loss=synergia.lattice.Diagnostics_loss("apertures_loss_%d.h5" % i,"aperture") 
               diag_loss.set_bunch(bunch_train_simulator.get_bunch_train().get_bunches()[i]) 
               stepper.get_lattice_simulator().get_lattice().add_loss_diagnostics(diag_loss)


# ofset bunches and adjust means 
      if bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm().has_this_rank():
            bunch_means=Core_diagnostics.calculate_mean(bunch_train_simulator.get_bunch_train().get_bunches()[i])
            particles = bunch_train_simulator.get_bunch_train().get_bunches()[i].get_local_particles()
            particles[:,0] = particles[:,0]-bunch_means[0]+opts.x_offset
            particles[:,1] = particles[:,1]-bunch_means[1]
            particles[:,2] = particles[:,2]-bunch_means[2]+opts.y_offset
            particles[:,3] = particles[:,3]-bunch_means[3]
            particles[:,4] = particles[:,4]-bunch_means[4]+opts.z_offset
            particles[:,5] = particles[:,5]-bunch_means[5]

            bunch_means=Core_diagnostics.calculate_mean(bunch_train_simulator.get_bunch_train().get_bunches()[i])
            bunch_stds=Core_diagnostics.calculate_std(bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means)
            bunch_mom2=Core_diagnostics.calculate_mom2(bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means)
            if bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm().get_rank()==0:
                print
                print "bunch # ",i," is perriodic=",bunch_train_simulator.get_bunch_train().get_bunches()[i].is_periodic() 
                print "bunch # ",i," length=",bunch_train_simulator.get_bunch_train().get_bunches()[i].get_z_period_length() 
                print "bunch # ",i," has longitudinal aperture=",bunch_train_simulator.get_bunch_train().get_bunches()[i].has_longitudinal_aperture() 
                print "bunch # ",i," number of real  particles= ",bunch_train_simulator.get_bunch_train().get_bunches()[i].get_real_num() 
                print "bunch # ",i," number of macroparticles= ",bunch_train_simulator.get_bunch_train().get_bunches()[i].get_total_num() 
                print "bunch # ",i," bucket index= ",bunch_train_simulator.get_bunch_train().get_bunches()[i].get_bucket_index() 
                print "bunch # ",i," initial offsets (x,xp,y,yp,ct,dpp)=(",bunch_means[0],", ",bunch_means[1],", ",bunch_means[2],", "\
                                        ,bunch_means[3],", ",bunch_means[4],", ",bunch_means[5] ,") [meters]" 
                print "bunch # ",i," initial stds (xrms,yrms,zrms)=(",bunch_stds[0],", ",bunch_stds[2],", ",beta*bunch_stds[4],") [meters]"                    
                print "bunch # ",i," :" 
                Core_diagnostics.print_bunch_parameters(bunch_mom2, beta)
                print "___________________________________________________________"          

#make propagator and propagete
propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(opts.checkpointperiod)
propagator.set_concurrent_io(opts.concurrentio)
propagator.propagate(bunch_train_simulator, opts.num_turns, opts.maxturns, opts.verbosity)           