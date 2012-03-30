#!/usr/bin/env python
import sys
import synergia
from booster_options import opts
import numpy as np
import string
from synergia.optics.one_turn_map import linear_one_turn_map
from mpi4py import MPI




#lattice=synergia.lattice.Mad8_reader().get_lattice("booster", "booster_0.5.lat")




#lattice_length = lattice.get_length()
#reference_particle= lattice.get_reference_particle()
#energy = reference_particle.get_total_energy()
#beta = reference_particle.get_beta()
#gamma = reference_particle.get_gamma()

#for elem in lattice.get_elements():
    #if elem.get_name() == "rf":
        #elem.set_double_attribute("volt", opts.rf_voltage)
        #elem.set_double_attribute("freq",  opts.freq/2.0*np.pi)
        #elem.set_double_attribute("lag",  0.)
        #print "type=",  elem.print_()
       ## elem=synergia.lattice.Lattice_element("rfcavity","rfb")


#for elem in lattice.get_elements():
    #if elem.get_type() == "multipole":
        #print "multipole is ", elem.get_type(), " atributes are", elem.get_double_attributes(), "string aributes are ", elem.get_string_attributes()


#if MPI.COMM_WORLD.Get_rank() ==0:
    #print "lattice name=",lattice.get_name()
    #print "lattice length=",lattice.get_length()
    #print "reference particle energy: ", energy
    #print "reference particle beta: ", beta
    #print "reference particle gamma: ", gamma
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
        if elem.get_type() == "rfcavity":
            elem.set_double_attribute("l", 0.)
            elem.set_double_attribute("volt", opts.rf_voltage)
            elem.set_double_attribute("freq",  freq)
            elem.set_double_attribute("lag",  0.)
            #print "type=",  elem.print_()

        name= elem.get_name()[0:4]
        if (name=="fmag"):
            elem.set_string_attribute("aperture_type","rectangular")
            elem.set_double_attribute("rectangular_aperture_width", 0.2)
            elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
        elif (name=="dmag"):
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
        if opts.set_corrections:
            if correction_dict.has_key(elem.get_name()):
               elem.set_double_attribute("k1", correction_dict[elem.get_name()])


        if elem.get_name() == "orf":
            #print "elem. type=",elem.get_type()
            elem.set_double_attribute("volt", opts.rf_voltage)
            elem.set_double_attribute("freq",  freq)
            elem.set_double_attribute("lag",  0.)


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
        if opts.set_corrections:
            if correction_dict.has_key(elem.get_name()):
               elem.set_double_attribute("k1", correction_dict[elem.get_name()])
              # print "element=",elem.print_()

        if elem.get_name() == "rf":
           # print "elem. type=",elem.get_type()
            elem.set_double_attribute("volt", opts.rf_voltage)
            elem.set_double_attribute("freq",  freq)
            elem.set_double_attribute("lag",  0.)


        name= elem.get_name()[0:4]
        if (name=="fmag"):
            #print "name=",elem.get_name()
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

else:
    raise RuntimeError,  "booster.py: unknown lattice file for booster "
    sys.exit(1)



lattice_simulator = synergia.simulation.Lattice_simulator(lattice, opts.map_order)
bunch_sp= lattice_simulator.get_bucket_length()
map = linear_one_turn_map(lattice_simulator)




energy = lattice.get_reference_particle().get_total_energy()
gamma = lattice.get_reference_particle().get_gamma()
if MPI.COMM_WORLD.Get_rank() ==0:
    print "lattice name=",lattice.get_name()
    print "lattice length=",lattice.get_length()
    print "reference particle energy= ", energy
    print "reference particle beta= ", beta
    print "reference particle gamma= ", gamma
    print "momentum p= ",lattice_simulator.get_lattice().get_reference_particle().get_four_momentum().get_momentum()
    print "freq=",freq

if MPI.COMM_WORLD.Get_rank() ==0:
    print " bucket_lenght=",bunch_sp
    print "lattice size=",len(lattice.get_elements())

[[ax,ay],[bx,by]] = synergia.optics.get_alpha_beta(map)
(tune_x, tune_y, tune_z) = synergia.optics.get_tunes(map)


if MPI.COMM_WORLD.Get_rank() ==0:
    print "Lattice functions assuming uncoupled map:"
    print "alpha x: ", ax
    print "alpha y: ", ay
    print "beta x: ", bx
    print "beta y: ", by
    print "tune x: ", tune_x, "(",1-tune_x,")"
    print "tune y: ", tune_y, "(",1-tune_y,")"
    print "tune z: ", tune_z, "(",1-tune_z,")"




tt=True
if  not((opts.latticefile=="booster_petrenko.lat") or ((opts.latticefile=="booster_Leo_sbend.lat") and opts.set_corrections )):
    elem= lattice.get_elements()[len(lattice.get_elements())-1]
    (betax, alphax, betay, alphay)=    (lattice_simulator.get_lattice_functions(elem).beta_x,\
    lattice_simulator.get_lattice_functions(elem).alpha_x,lattice_simulator.get_lattice_functions(elem).beta_y,\
    lattice_simulator.get_lattice_functions(elem).alpha_y)
    if MPI.COMM_WORLD.Get_rank() ==0:
       print "lattice functions for last element ",elem.get_name(), "(alpha_x, alpha_y, beta_x, beta_y) = (%g, %g, %g, %g)" % (alphax, alphay, betax, betay)
       print " nu_h=", lattice_simulator.get_horizontal_tune(), "p= ",lattice_simulator.get_lattice().get_reference_particle().get_four_momentum().get_momentum()
       print " nu_v=", lattice_simulator.get_vertical_tune()



# now make the collective operators

impedance=opts.impedance
if impedance:
    zgrid=40
    imped_f= synergia.collective.Impedance("LamF_wake.dat",lattice_length, bunch_sp,zgrid, "rectangular",60)
    if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "WAKES FOR F MAGNET read from ", imped_f.get_wake_file_name()
        print "F mag orbith length=",imped_f.get_orbit_length()
        print "F mag z_grid=",imped_f.get_z_grid()
        print "F mag stored turns=",imped_f.get_nstored_turns()

    imped_d= synergia.collective.Impedance("LamD_wake.dat",lattice_length, bunch_sp,zgrid, "rectangular",60)
    if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "WAKES FOR F MAGNET read from", imped_d.get_wake_file_name()
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

    #grid_shape_L=opts.scgrid[:]
    grid_shape_L=[128,128,64]
    radius_L=opts.get("apertureL")
    pipe_size_L=[2.*radius_L, 2.*radius_L, lattice_simulator.get_bucket_length()]
    spc_else=synergia.collective.Space_charge_rectangular(pipe_size_L, grid_shape_L)
    if MPI.COMM_WORLD.Get_rank() ==0:
        print
        print "pipe_size L section=",spc_else.get_pipe_size()
        print "grid for spch L section=", spc_else.get_grid_shape()
        print "___________________________________________________________"


# now make stepper
if impedance or space_charge :
    list_choice_map=synergia.simulation.List_choice_map()
    operators_fmag=[]
    steps_per_fmag=1
    if impedance: operators_fmag.append(imped_f)
    if space_charge: operators_fmag.append(spc_f)
    kicks_fmag=synergia.simulation.Kicks(operators_fmag,steps_per_fmag)

    operators_dmag=[]
    steps_per_dmag=1
    if impedance: operators_dmag.append(imped_d)
    if space_charge: operators_dmag.append(spc_d)
    kicks_dmag=synergia.simulation.Kicks(operators_dmag,steps_per_dmag)
    for elem in lattice.get_elements():
        if  (opts.latticefile=="booster_0.5.lat") or  (opts.latticefile=="booster_classic.lat"):
            name= elem.get_name()[0:4]
        elif  (opts.latticefile=="booster_petrenko.lat") or (opts.latticefile=="booster_petrenko_zero.lat") or \
                 (opts.latticefile=="booster_petrenko_leo.lat"):
              name= elem.get_name()[1:5]
        else:
            raise RuntimeError,  "booster.py: before list_choice_map: unknown lattice file for booster "
            sys.exit(1)

        if (name=="fmag"): list_choice_map[elem.get_name()]=kicks_fmag
        if (name=="dmag"): list_choice_map[elem.get_name()]=kicks_dmag

    operators_else=[]
    num_steps_else=opts.num_steps_else
    if space_charge: operators_else.append(spc_else)
    kicks_else=synergia.simulation.Kicks(operators_else,num_steps_else)
    list_choice_map["else"]=kicks_else

    stepper=synergia.simulation.Split_operator_stepper_choice(lattice_simulator,list_choice_map, True)



else:
    no_op = synergia.simulation.Dummy_collective_operator("stub")
    stepper = synergia.simulation.Split_operator_stepper(
                            lattice_simulator, no_op, opts.num_steps)
#stepper.print_()



# now create bunches

rms_index=[0,2,4]
arms=opts.xrms #np.sqrt(opts.emit*beta_x)
brms=opts.yrms #np.sqrt(opts.emit*beta_y)
crms=opts.zrms
if MPI.COMM_WORLD.Get_rank() ==0:
    print "input xrms=",arms
    print "input yrms=",brms
    print "input zrms=",crms

#covar = synergia.optics.matching._get_correlation_matrix(map,arms,brms,crms,beta,rms_index)
num_bunches=opts.num_bunches
bunch_diag_train=synergia.bunch.Bunch_with_diagnostics_train(num_bunches,bunch_sp, synergia.utils.Commxx())
for bunchnum in range(0,num_bunches):
    if bunch_diag_train.is_on_this_rank(bunchnum):
        commxx=bunch_diag_train.get_comm(bunchnum)
        bunch= synergia.optics.generate_matched_bunch(lattice_simulator,
                                                arms,brms,crms,
                                                opts.num_real_particles,#*(bunchnum+1),
                                                opts.num_macro_particles,rms_index,
                                                seed=opts.seed, bunch_index=bunchnum,comm=commxx, periodic=True)
        particles = bunch.get_local_particles()
        # apply offset to bunch
        particles[:,0] = particles[:,0]+opts.x_offset
        particles[:,2] = particles[:,2]+opts.y_offset
        particles[:,4] = particles[:,4]+opts.z_offset
        diagnostics_actions = synergia.simulation.Diagnostics_actions()
        bunch_diag=synergia.bunch.Bunch_with_diagnostics(bunch, diagnostics_actions)
        bunch_diag.add_per_step_diagnostics(synergia.bunch.Diagnostics_full2(bunch, "step_full2-%02d.h5"%bunchnum))
        bunch_diag.add_per_turn_diagnostics(synergia.bunch.Diagnostics_particles(bunch,"turn_particles-%02d.h5"%bunchnum,0,0,100))
        bunch_diag_train.set_bunch_diag_sptr(bunchnum, bunch_diag)
        if space_charge:
            spc_f.set_fftw_helper(commxx)
            spc_d.set_fftw_helper(commxx)
            spc_else.set_fftw_helper(commxx)

        real_num=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_real_num()
        macro_num=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_total_num()
        bucket_index=bunch_diag_train.get_bunch_diag_sptr(bunchnum).get_bunch_sptr().get_bucket_index()
        if commxx.get_rank() ==0:
            print "bunch # ",bunchnum ," number of real  particles= ",real_num
            print "bunch # ",bunchnum ," number of macroparticles= ",macro_num
            print "bunch # ",bunchnum ," bucket index", bucket_index
            print "___________________________________________________________"
if MPI.COMM_WORLD.Get_rank() ==0:
    print "train bunch space=",bunch_diag_train.get_bunch_separation()






propagator = synergia.simulation.Propagator(stepper)
propagator.propagate(bunch_diag_train, opts.num_turns, opts.verbose)