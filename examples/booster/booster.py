#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import synergia
from booster_options import opts
import numpy as np
import string
from synergia.optics.one_turn_map import linear_one_turn_map
from synergia.utils import Commxx
from mpi4py import MPI



try:
  quad_correctors_H=[] 
  quad_correctors_V=[] 
  sextupole_correctors_H=[] 
  sextupole_correctors_V=[]  
  
  
  if  (opts.latticefile=="booster_0.5.lat") or  (opts.latticefile=="booster_classic.lat"):  
      lattice= synergia.lattice.Lattice("booster")   
      rfcells = [14,15,16,17,19,21,22,23,24]
      cell_line_no_rf = range(0,25)
      cell_line = range(0,25)
      for cell in range(1,25):
	# print "bcell=","bcel%02d"%cell
	  cell_line_no_rf[cell]=synergia.lattice.Mad8_reader().get_lattice("bcel%02d"%cell ,opts.latticefile)
	  cell_line[cell]=synergia.lattice.Lattice("bcel%02d"%cell) 
      
      #lattice_cel=synergia.lattice.Lattice("bcel") 
      icell=14
      
      drift1_length=0.22
      drift2_length=2.56
      for cell in range(1,25):
	  if cell in rfcells:
	      for elem in cell_line_no_rf[cell].get_elements():
		  
		  if opts.chef_propagate:
		      elem.set_string_attribute("extractor_type", "chef_propagate")
		  elif  opts.chef_map:   
		      elem.set_string_attribute("extractor_type", "chef_map")
		      
		  if elem.get_name() == "longa":
		      longa1=synergia.lattice.Lattice_element("drift", "longa_1")
		      longa1.set_double_attribute("l", drift1_length)
		      lattice.append(longa1)
		      cell_line[cell].append(longa1)
		      
		      rfa=synergia.lattice.Lattice_element("rfcavity", "rfa")
		      #rfa.set_double_attribute("volt", opts.rf_voltage)
		      #rfa.set_double_attribute("freq",  opts.freq)
		      #rfa.set_double_attribute("lag",  0.)
		      lattice.append(rfa)
		      cell_line[cell].append(rfa)
		      
		      
		      longa2=synergia.lattice.Lattice_element("drift", "longa_2")
		      longa2.set_double_attribute("l", drift2_length)
		      lattice.append(longa2)
		      cell_line[cell].append(longa2)
		      
		      rfb=synergia.lattice.Lattice_element("rfcavity", "rfb")
		      #rfb.set_double_attribute("volt", opts.rf_voltage)
		      #rfb.set_double_attribute("freq",  opts.freq)
		      #rfb.set_double_attribute("lag",  0.)
		      lattice.append(rfb)
		      cell_line[cell].append(rfb)
		      
		      
		      longa3=synergia.lattice.Lattice_element("drift", "longa_3")
		      longa3.set_double_attribute("l", drift1_length)
		      lattice.append(longa3)
		      cell_line[cell].append(longa3)
		  else:     
		      lattice.append(elem)
		      cell_line[cell].append(elem)
	  else:
	      for elem in cell_line_no_rf[cell].get_elements():
		  lattice.append(elem)  
		  cell_line[cell].append(elem) 
			
	  cell_line[cell].set_reference_particle(cell_line_no_rf[1].get_reference_particle())
      
      lattice.set_reference_particle(cell_line_no_rf[1].get_reference_particle())
      lattice_length=lattice.get_length() 
      reference_particle = lattice.get_reference_particle()
      beta = reference_particle.get_beta()
    
      
      harmon=opts.num_buckets
      freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
      #freq=37867099.7584
      
      print "freq=",freq
      for elem in lattice.get_elements():
	  if opts.chef_propagate:
	      elem.set_string_attribute("extractor_type", "chef_propagate")
	
	  if elem.get_type() == "rfcavity":
	      elem.set_double_attribute("l", 0.)
	      elem.set_double_attribute("volt", opts.rf_voltage)
	      elem.set_double_attribute("freq",  freq)
	      elem.set_double_attribute("lag",  0.)
	    # print "type=",  elem.print_()         
	  if opts.aperture:
	      name= elem.get_name()[0:4]
	      if (name=="fmag"):
		# print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
	      elif (name=="dmag"):
		  #print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureD"))
	      else:
		  elem.set_string_attribute("aperture_type","circular")
		  elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL"))
		    
  elif (opts.latticefile=="booster_petrenko.lat") or (opts.latticefile=="booster_petrenko_zero.lat"):
      lattice=synergia.lattice.Mad8_reader().get_lattice("machine", opts.latticefile)
      
      lattice_length=lattice.get_length() 
      reference_particle = lattice.get_reference_particle()
      beta = lattice.get_reference_particle().get_beta()
      harmon=opts.num_buckets
      freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
      #freq=37867099.7584
      
      if opts.set_corrections:
	correction_dict={} 
	cfile = open(opts.correction_file,"r")
	for read_line in cfile.readlines():
	      split_line=read_line.split()
	      if ( (len(split_line)==6) and (split_line[3]=="QUAD") ) :
		  correction_dict[split_line[0].swapcase()]=eval(split_line[2])
		# print "dict =",correction_dict
      

      for elem in lattice.get_elements():
	  if opts.chef_propagate:
	      elem.set_string_attribute("extractor_type", "chef_propagate")
	  elif  opts.chef_map:   
	      elem.set_string_attribute("extractor_type", "chef_map")
	      
	  if opts.set_corrections:
	      if correction_dict.has_key(elem.get_name()):
		elem.set_double_attribute("k1", correction_dict[elem.get_name()])
	      
	      
	  if elem.get_name() == "orf":
	      #print "elem. type=",elem.get_type()
	      elem.set_double_attribute("volt", opts.rf_voltage)
	      elem.set_double_attribute("freq",  freq)
	      elem.set_double_attribute("lag",  0.) 
	    
	  if opts.aperture:   
	      name= elem.get_name()[1:5]
	      if (name=="fmag"):
		# print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
	      elif (name=="dmag"):
		# print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureD"))
	      else:
		  elem.set_string_attribute("aperture_type","circular")
		  elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL"))   
	      
	  #if elem.get_name() =="kpinger":
	      #print " UUUUUU name=",elem.get_name(),"  type=",elem.get_type()
	      #print "element=",elem.print_() 
  elif (opts.latticefile=="booster_2012.lat") or (opts.latticefile=="booster_newtunes.lat"): # or (opts.latticefile=="booster_expchrom.lat"):"  
      lattice=synergia.lattice.Mad8_reader().get_lattice("machine", opts.latticefile)

	    
      lattice_length=lattice.get_length() 
      reference_particle = lattice.get_reference_particle()
      beta = lattice.get_reference_particle().get_beta()
      harmon=opts.num_buckets
      freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
      #freq=37867099.7584
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print "initial frequency=",freq      

      for elem in lattice.get_elements():
	  if opts.chef_propagate:
	      elem.set_string_attribute("extractor_type", "chef_propagate")
	  elif  opts.chef_map:   
	      elem.set_string_attribute("extractor_type", "chef_map")
				
	      
	  if elem.get_name() == "arf":
	      elem.set_double_attribute("volt", opts.rf_voltage)
	      elem.set_double_attribute("freq",  freq)
	      elem.set_double_attribute("lag",  0.) 
	    # print "element=",elem.print_()
	      
	  if opts.aperture:   
	      name= elem.get_name()[0:5]
	      if (name=="bfmag"):
		# print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
	      elif (name=="bdmag"):
		# print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureD"))
	      else:
		  elem.set_string_attribute("aperture_type","circular")
		  elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL"))   
	  
	  names= elem.get_name()[0:4] 
	  if ((names=="ssxs")  and (elem.get_type()=="sextupole")):        
	    sextupole_correctors_H.append(elem) 
	    #print "name=",elem.get_name() 
	  if ((names=="ssxl")   and (elem.get_type()=="sextupole")):        
	    sextupole_correctors_V.append(elem) 
	    #print "name=",elem.get_name()
	  nameq= elem.get_name()[0:3] 
	  nameql= elem.get_name()[0:6]
	  if ((nameq=="qql")  and (elem.get_type()=="quadrupole") and (nameql != "qqlerr") ):                    
	    quad_correctors_H.append(elem)     
	
	  if ((nameq=="qqs")  and (elem.get_type()=="quadrupole") and (nameql != "qqserr")):                    
	    quad_correctors_V.append(elem)        
	  #if elem.get_name() =="khpinger":
	      #print " UUUUUU name=",elem.get_name(),"  type=",elem.get_type()
	      #print "element=",elem.print_()   
	
  elif (opts.latticefile=="booster_Leo_sbend.lat"):    
      lattice=synergia.lattice.Mad8_reader().get_lattice("booster", opts.latticefile)
      lattice_length=lattice.get_length() 
      reference_particle = lattice.get_reference_particle()
      beta = reference_particle.get_beta()
      harmon=opts.num_buckets
      freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
      
      if opts.set_corrections:
	correction_dict={} 
	cfile = open(opts.correction_file,"r")
	for read_line in cfile.readlines():
	      split_line=read_line.split()
	      if ( (len(split_line)==6) and (split_line[3]=="QUAD") ) :
		  correction_dict[split_line[0].swapcase()[1:]]=eval(split_line[2])
		# print "dict =",correction_dict
		  
		  
	
      for elem in lattice.get_elements():
	  if opts.chef_propagate:
	      elem.set_string_attribute("extractor_type", "chef_propagate")
	  elif  opts.chef_map:   
	      elem.set_string_attribute("extractor_type", "chef_map")
	  
	  if opts.set_corrections:
	      if correction_dict.has_key(elem.get_name()):
		elem.set_double_attribute("k1", correction_dict[elem.get_name()])
		# print "element=",elem.print_() 
	      
	  if (elem.get_name() == "rf") :            
	    # print "elem. type=",elem.get_type()
	      elem.set_double_attribute("volt", opts.rf_voltage)
	      elem.set_double_attribute("freq",  freq)
	      elem.set_double_attribute("lag",  0.) 
	    
	    
	  if opts.aperture:   
	      name= elem.get_name()[0:4]
	      if (name=="fmag"):                
		  #print "name=",elem.get_name()                
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
		
	      elif (name=="dmag"):
		  #print "name=",elem.get_name()
		  elem.set_string_attribute("aperture_type","rectangular")
		  elem.set_double_attribute("rectangular_aperture_width", 0.2)
		  elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureD"))
		  
		
	      else:
		  elem.set_string_attribute("aperture_type","circular")
		  elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL")) 
				  
	      #elem.set_string_attribute("aperture_type","circular")
	      #elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL"))    
		  
  else:
      print "latticefile=", opts.latticefile
      raise RuntimeError,  "booster.py: unknown lattice file for booster " 
      sys.exit(1)   

    
  #for qua in quad_correctors_H:
      #print " quad_correctors_H=",  qua.get_name() 
  #print 
  #for qua in quad_correctors_V:
      #print " quad_correctors_V=",  qua.get_name() 
	    
	    
  lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)  
  lattice_simulator.get_lattice().get_reference_particle()        
  bunch_sp= lattice_simulator.get_bucket_length()

  if opts.adjust_tunes:
    lattice_simulator.adjust_tunes(opts.tuneH, opts.tuneV ,quad_correctors_H, quad_correctors_V)

  
  if opts.adjust_chromaticity:
      lattice_simulator.adjust_chromaticities(opts.chromH, opts.chromV, \
	  sextupole_correctors_H,sextupole_correctors_V)




  

    
  # now make the collective operators

  impedance=opts.impedance
  if impedance:
      zgrid=1000
      wn=opts.wave_number[:]
      imped_f=synergia.collective.Impedance(opts.wakefileF, "XLXTYLYTZpp", zgrid, lattice_length, bunch_sp,opts.registred_turns, opts.full_machine,wn); 
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print
	  print "WAKES FOR F MAGNET read from ", imped_f.get_wake_field().get_wake_file_name()
	  print "F mag orbith length=",imped_f.get_orbit_length()
	  print "F mag z_grid=",imped_f.get_z_grid()
	  print "F mag stored turns=",imped_f.get_nstored_turns()
      imped_d=synergia.collective.Impedance(opts.wakefileD, "XLXTYLYTZpp", zgrid, lattice_length, bunch_sp,opts.registred_turns, opts.full_machine,wn)    
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print
	  print "WAKES FOR D MAGNET read from", imped_d.get_wake_field().get_wake_file_name()
	  print "D mag orbith length=",imped_d.get_orbit_length()
	  print "D mag z_grid=",imped_d.get_z_grid()
	  print "D mag stored turns=",imped_d.get_nstored_turns()
	  print "___________________________________________________________"
  space_charge=opts.space_charge
  if space_charge:
      grid_shape_f=opts.scgrid[:]
      radiusx_f=0.1
      radiusy_f= opts.get("apertureF")  
      grid_shape_f[0] *= int(radiusx_f/opts.grid_ref_distance)
      grid_shape_f[1] *= int(radiusy_f/opts.grid_ref_distance)
      pipe_size_f=[2.*radiusx_f, 2.*radiusy_f, lattice_simulator.get_bucket_length()]
      spc_f=synergia.collective.Space_charge_rectangular(pipe_size_f, grid_shape_f)
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print
	  print "pipe_size F magnet=",spc_f.get_pipe_size()
	  print "grid for spch F magnet=",spc_f.get_grid_shape() 
	
  #*********************************

      grid_shape_d=opts.scgrid[:]
      radiusx_d=0.1
      radiusy_d= opts.get("apertureD") 
      
      grid_shape_d[0] *= int(radiusx_d/opts.grid_ref_distance)
      grid_shape_d[1] *= int(radiusy_d/opts.grid_ref_distance)
      pipe_size_d=[2.*radiusx_d, 2.*radiusy_d, lattice_simulator.get_bucket_length()]
      spc_d=synergia.collective.Space_charge_rectangular(pipe_size_d, grid_shape_d)
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print
	  print "pipe_size D magnet=",spc_d.get_pipe_size()
	  print "grid for spch D magnet=",spc_d.get_grid_shape() 
	
      grid_shape_L=opts.scgrid_L[:]
      radius_L=opts.get("apertureL")     
      pipe_size_L=[2.*radius_L, 2.*radius_L, lattice_simulator.get_bucket_length()]        
      spc_else=synergia.collective.Space_charge_rectangular(pipe_size_L, grid_shape_L) 
      if MPI.COMM_WORLD.Get_rank() ==0:
	  print
	  print "pipe_size L section=",spc_else.get_pipe_size()
	  print "grid for spch L section=", spc_else.get_grid_shape() 
	  print "___________________________________________________________"

	  
	  
  # now make stepper  
  bpms=opts.bpms
  if bpms:
      bpm_measure = synergia.simulation.Dummy_collective_operator("bmp_measure")
  if impedance or space_charge or bpms :    
      list_choice_map=synergia.simulation.List_choice_map()
      operators_fmag=[]
      steps_per_fmag=opts.steps_per_fmag
      if impedance: operators_fmag.append(imped_f)
      if space_charge: operators_fmag.append(spc_f)
      kicks_fmag=synergia.simulation.Kicks(operators_fmag,steps_per_fmag)
      
      operators_dmag=[]
      steps_per_dmag=opts.steps_per_dmag    
      if impedance: operators_dmag.append(imped_d)
      if space_charge: operators_dmag.append(spc_d)
      kicks_dmag=synergia.simulation.Kicks(operators_dmag,steps_per_dmag)
      
      steps_per_bpm=1
      operators_bpm=[]
      if bpms: operators_bpm.append(bpm_measure)
      kicks_bpm=synergia.simulation.Kicks(operators_bpm,steps_per_bpm)
      
      for elem in lattice.get_elements():
	  if  (opts.latticefile=="booster_0.5.lat") or  (opts.latticefile=="booster_classic.lat")or \
		  (opts.latticefile=="booster_Leo_sbend.lat"):
	      name= elem.get_name()[0:4]
	      name1=elem.get_name()[0:3]
	  elif  (opts.latticefile=="booster_petrenko.lat") or (opts.latticefile=="booster_petrenko_zero.lat") or \
		  (opts.latticefile=="booster_petrenko_leo.lat"):
		name= elem.get_name()[1:5]
		name1=elem.get_name()[0:3]
	  elif  (opts.latticefile=="booster_2012.lat") or (opts.latticefile=="booster_newtunes.lat"): 
		name= elem.get_name()  
		name1= elem.get_name()
	  else:
	      raise RuntimeError,  "booster.py: before list_choice_map: unknown lattice file for booster " 
	      sys.exit(1)
	  if impedance or space_charge:  
	      if  (opts.latticefile=="booster_2012.lat") or (opts.latticefile=="booster_zerochrom.lat") or \
		  (opts.latticefile=="booster_expchrom.lat"):
		  if (name=="bfmag"): list_choice_map[elem.get_name()]=kicks_fmag
		  if (name=="bdmag"): list_choice_map[elem.get_name()]=kicks_dmag
	      else:                            
		if (name=="fmag"): list_choice_map[elem.get_name()]=kicks_fmag
		if (name=="dmag"): list_choice_map[elem.get_name()]=kicks_dmag
	  if bpms:
	      if  (opts.latticefile=="booster_2012.lat") or (opts.latticefile=="booster_zerochrom.lat") or \
		  (opts.latticefile=="booster_expchrom.lat"):
		  if (name=="obpm"): 
		      list_choice_map[elem.get_name()]=kicks_bpm
		      if MPI.COMM_WORLD.Get_rank() ==0:
			  print "BPMS are:  ",elem.get_name()                  
	      else:                             
		  if (name1=="bpm"): 
		      list_choice_map[elem.get_name()]=kicks_bpm
		      if MPI.COMM_WORLD.Get_rank() ==0:
			  print "BPMS are:  ",elem.get_name()
      
      operators_else=[]
      num_steps_else=opts.num_steps_else 
      #no_op = synergia.simulation.Dummy_collective_operator("stub")
      #if space_charge: operators_else.append(no_op) 
      if space_charge: operators_else.append(spc_else)    
      kicks_else=synergia.simulation.Kicks(operators_else,num_steps_else)
      list_choice_map["else"]=kicks_else
      if impedance or space_charge:
	  if MPI.COMM_WORLD.Get_rank() ==0:
	      print
	      print " steps_per_fmag=",steps_per_fmag
	      print " steps_per_dmag=",steps_per_dmag
	      print " num_steps_else=",num_steps_else
	      print "___________________________________________________________"
      
      stepper=synergia.simulation.Split_operator_stepper_choice(lattice_simulator,list_choice_map, True)
  else:   
      stepper = synergia.simulation.Independent_stepper(lattice_simulator, opts.num_steps)
      if MPI.COMM_WORLD.Get_rank() ==0:                                                   
	  print   "no collective effects, no bpm, independent operator"                     
	  #stepper.print_()  
     
      #stepper = synergia.simulation. Independent_stepper_elements(lattice_simulator, 1)
      #if MPI.COMM_WORLD.Get_rank() ==0:                                                   
	  #print   "independenet stepper elements"                     
	  #stepper.print_()            
	  
    
      #no_op = synergia.simulation.Dummy_collective_operator("stub")
      #stepper = synergia.simulation.Split_operator_stepper(
			      #lattice_simulator, no_op, opts.num_steps) 
      #if MPI.COMM_WORLD.Get_rank() ==0:                                                   
	  #print   "no collective effects, no bpm, split_operator_stepper"                     
	  #stepper.print_()                    
  #if MPI.COMM_WORLD.Get_rank() ==0:
     #stepper.print_()
     
    
  
	   
  energy = lattice.get_reference_particle().get_total_energy()
  gamma = lattice.get_reference_particle().get_gamma()
  if MPI.COMM_WORLD.Get_rank() ==0:
      print "lattice name=",lattice.get_name() 
      print "lattice length=",lattice.get_length()  
      print "reference particle energy= ", energy
      print "reference particle beta= ", beta
      print "reference particle gamma= ", gamma
      print "momentum p= ",stepper.get_lattice_simulator().get_lattice().get_reference_particle().get_four_momentum().get_momentum()
      print "freq=",freq

  if MPI.COMM_WORLD.Get_rank() ==0:
      print " bucket_lenght=",bunch_sp
      print "lattice size=",len(lattice.get_elements())

  #mapold = linear_one_turn_map(stepper.get_lattice_simulator())
  #print "map old=",mapold
  map= stepper.get_lattice_simulator().get_linear_one_turn_map()
  #print "map lattice sim=",map
  [[ax,ay],[bx,by]] = synergia.optics.get_alpha_beta(map)
  (tune_x, tune_y, tune_z) = synergia.optics.get_tunes(map)


  #stepper.get_lattice_simulator().print_lattice_functions()
  #chef_frac_tunex=stepper.get_lattice_simulator().get_horizontal_tune()
  #chef_frac_tuney=stepper.get_lattice_simulator().get_vertical_tune()
  #chef_eigen_tunex=stepper.get_lattice_simulator().get_horizontal_tune(True)
  #chef_eigen_tuney=stepper.get_lattice_simulator().get_vertical_tune(True)
  #horizontal_chromaticity=stepper.get_lattice_simulator().get_horizontal_chromaticity()
  #vertical_chromaticity=stepper.get_lattice_simulator().get_vertical_chromaticity()
  #momentum_compaction=stepper.get_lattice_simulator().get_momentum_compaction()
  #slip_factor=stepper.get_lattice_simulator().get_slip_factor()
  if MPI.COMM_WORLD.Get_rank() ==0:
      print "Lattice functions assuming uncoupled map:"
      print "alpha x: ", ax
      print "alpha y: ", ay
      print "beta x: ", bx
      print "beta y: ", by
      #print "chef FracTune x: ", chef_frac_tunex, ", EigenTune x: ", chef_eigen_tunex, ", map tune x: ", tune_x, "(",1-tune_x,")"
      #print "chef FracTune y: ", chef_frac_tuney, ", EigenTune y: ", chef_eigen_tuney, ", map tune y: ", tune_y, "(",1-tune_y,")"
      #print "                           map tune z: ", tune_z, "(",1-tune_z,")"  
      #print " horizontal chromaticity: ", horizontal_chromaticity
      #print " vertical   chromaticity: ", vertical_chromaticity
      #print " momentum compaction: ", momentum_compaction
      #print " slip factor: ", slip_factor





  # now create bunches
    
  rms_index=[0,2,4]
  arms=opts.xrms #np.sqrt(opts.emit*beta_x)
  brms=opts.yrms #np.sqrt(opts.emit*beta_y)
  crms=opts.zrms  # z=beta*c*t!
  if MPI.COMM_WORLD.Get_rank() ==0:
      print "input xrms=",arms
      print "input yrms=",brms
      print "input zrms=",crms


  num_bunches=opts.num_bunches
  bunches = []
  parent_comm = Commxx()
  comms = synergia.utils.generate_subcomms(parent_comm, opts.num_bunches)
  for i in range(0, num_bunches):
      commx=comms[i]
      bunch= synergia.optics.generate_matched_bunch(stepper.get_lattice_simulator(),
						  arms,brms,crms,
						  opts.num_real_particles,
						  opts.num_macroparticles,rms_index,
						  seed=opts.seed, bunch_index=i,comm=commx, periodic=opts.periodic)
      particles = bunch.get_local_particles()
      particles[:,0] = particles[:,0]+opts.x_offset*np.cos(2.*np.pi*i/float(num_bunches))
      particles[:,2] = particles[:,2]+opts.y_offset*np.cos(2.*np.pi*i/float(num_bunches))
      particles[:,4] = particles[:,4]+opts.z_offset*np.cos(2.*np.pi*i/float(num_bunches))
      bunches.append(bunch)
     


  bunch_train = synergia.bunch.Bunch_train(bunches, bunch_sp)
  bunch_train_simulator = synergia.simulation.Bunch_train_simulator(bunch_train)


  Step_Numbers=[] 
  step_number=0
  for step  in stepper.get_steps():
    step_number += 1
    for op in step.get_operators():	
	if(op.get_name()=="bmp_measure"):
	    Step_Numbers.append(step_number)

  for i in range(0, bunch_train.get_size()):        
      bunch_train_simulator.add_per_step(i, synergia.bunch.Diagnostics_full2("step_full2_%d.h5" % i), Step_Numbers)  
      bunch_train_simulator.add_per_turn(i, synergia.bunch.Diagnostics_particles("turn_particles_%d.h5" % i), opts.turn_period)

      if space_charge:
        spc_optim=opts.spc_optim
	if bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm().has_this_rank():
	  #comm_spc=bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm()
	  comm_spc= synergia.utils.make_optimal_spc_comm(bunch_train_simulator.get_bunch_train().get_bunches()[i].get_comm(),spc_optim);
	  spc_f.set_fftw_helper(comm_spc)
	  spc_d.set_fftw_helper(comm_spc)
	  spc_else.set_fftw_helper(comm_spc)
	 

      real_num=bunch_train_simulator.get_bunch_train().get_bunches()[i].get_real_num()
      macro_num= bunch_train_simulator.get_bunch_train().get_bunches()[i].get_total_num()
      bucket_index= bunch_train_simulator.get_bunch_train().get_bunches()[i].get_bucket_index()
      #print "bunch # ",i ," num_local_particles=",bunch_train_simulator.get_bunch_train().get_bunches()[i].get_local_num(),\
               #"  rank= ",MPI.COMM_WORLD.Get_rank()
      if  MPI.COMM_WORLD.Get_rank() ==0:
	      print "bunch # ",i ," number of real  particles= ",real_num
	      print "bunch # ",i ," number of macroparticles= ",macro_num
	      print "bunch # ",i ," bucket index", bucket_index          
	      print "___________________________________________________________"
  if num_bunches >1:
    if MPI.COMM_WORLD.Get_rank() ==0:
      print "train bunch space=",bunch_train_simulator.get_bunch_train().get_spacings()[0]
      print  
       




  #bunch_diag_train=synergia.bunch.Bunch_with_diagnostics_train(num_bunches,bunch_sp, MPI.COMM_WORLD)
  #for bunchnum in range(0,num_bunches):
      #if bunch_diag_train.is_on_this_rank(bunchnum):
	  #commx=bunch_diag_train.get_comm(bunchnum)
	  ##bunch = synergia.optics.generate_matched_bunch_normalcoo(stepper.get_lattice_simulator(),
						  ##actions,
						  ##opts.num_real_particles,
						  ##opts.num_macro_particles,
						  ##seed=opts.seed, bunch_index=bunchnum,comm=commx, periodic=opts.periodic)
		  
	  #bunch= synergia.optics.generate_matched_bunch(stepper.get_lattice_simulator(),
						  #arms,brms,crms,
						  #opts.num_real_particles,
						  #opts.num_macro_particles,rms_index,
						  #seed=opts.seed, bunch_index=bunchnum,comm=commx, periodic=opts.periodic)
	  #particles = bunch.get_local_particles()
	  ## apply offset to bunch
	  #particles[:,0] = particles[:,0]+opts.x_offset
	  #particles[:,2] = particles[:,2]+opts.y_offset
	  #particles[:,4] = particles[:,4]+opts.z_offset 
	  #diagnostics_actions = synergia.simulation.Standard_diagnostics_actions()
	  #bunch_diag=synergia.bunch.Bunch_with_diagnostics(bunch, diagnostics_actions)
	  #bunch_diag.add_per_step_diagnostics(synergia.bunch.Diagnostics_full2(bunch, "step_full2-%02d.h5"%bunchnum))        
	  #bunch_diag.add_per_turn_diagnostics(synergia.bunch.Diagnostics_particles(bunch,"turn_particles-%02d.h5"%bunchnum,0,0,opts.turn_write))
	  #if bpms:  
	      #diagnostics_actions.addto_list_with_steps_for_diagnostics("bpm")
	  #else:
	      #diagnostics_actions.addto_list_with_steps_for_diagnostics("all_steps")
	  #bunch_diag_train.set_bunch_diag_sptr(bunchnum, bunch_diag)
	  #if space_charge:
	      #spc_f.set_fftw_helper(commx)
	      #spc_d.set_fftw_helper(commx)
	      #spc_else.set_fftw_helper(commx)
	      
	  #real_num=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_real_num()
	  #macro_num=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_total_num()
	  #bucket_index=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_bucket_index()
	##  if ((bucket_index >= ac_bucket_min) and (bucket_index <= ac_bucket_max)): action_on_bunch[bunchnum]=1    
	  #if commx.Get_rank() ==0:
	      #print "bunch # ",bunchnum ," number of real  particles= ",real_num
	      #print "bunch # ",bunchnum ," number of macroparticles= ",macro_num
	      #print "bunch # ",bunchnum ," bucket index", bucket_index
	      ##print "action_on_bunch=", action_on_bunch
	      #print "___________________________________________________________"
  #if MPI.COMM_WORLD.Get_rank() ==0:
      #print "train bunch space=",bunch_diag_train.get_bunch_separation()


  
  propagator = synergia.simulation.Propagator(stepper)
  propagator.set_checkpoint_period(opts.checkpointperiod)
  propagator.set_concurrent_io(opts.concurrentio)
  propagator.propagate(bunch_train_simulator, opts.num_turns, opts.maxturns, opts.verbosity)
except Exception, e:
    sys.stderr.write(str(e) + '\n')
    MPI.COMM_WORLD.Abort(777)
