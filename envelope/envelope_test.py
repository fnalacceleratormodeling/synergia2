import os.path
import sys
import synergia



if ( __name__ == '__main__'):
  
    myopts = synergia.Options("envelope_match")
#    myopts.add("latticefile","fodo.lat","",str)
    myopts.add("latticefile","fobodobo_s.lat","",str)
    myopts.add("emitx",8.89533703303356e-7,"X emittance",float)
    myopts.add("emity",1.89533703303356e-5,"Y emittance",float)
    myopts.add("energy",9.004401675138,"",float)
    myopts.add("bunchnp",5.0e+11,"number of particles per bunch",float)
    
    
    charge = 1.0  # electron charge in C
    initial_phase = 0.0
    scaling_frequency = 47713451.5923694
    
    energy = myopts.get("energy")
    mass = synergia.PH_NORM_mp
    kinetic_energy = energy-mass
    emittance = myopts.get("emittance")
    emitx=myopts.get("emitx")
    emity=myopts.get("emity")
    
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,"model",kinetic_energy,
                        scaling_frequency)
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)	
    beta = beam_parameters.get_beta()			 
 
    line_length = gourmet.orbit_length()
    bunch_spacing = line_length/588.0		
    current = myopts.get("bunchnp")* \
        synergia.physics_constants.PH_MKS_e/ \
        (bunch_spacing/(beta*synergia.physics_constants.PH_MKS_c))

		
     
    #[sigma_x, sigma_xprime, r_x, sigma_y,  sigma_yprime,  r_y] =synergia.matching.envelope_match_alex(emitx,emity,current,gourmet)
    #print " "
    #print " after matching, printed from  main"
    #print " "
    #print "sigma_x = %5.8f" % sigma_x
    #print "sigma_xprime= %5.8f" % sigma_xprime
    #print "r_x= %5.8f" % r_x
    #print "sigma_y = %5.8f" % sigma_y
    #print "sigma_yprime= %5.8f " % sigma_yprime
    #print "r_y= %5.8f" % r_y
    
    
   
    
    retval=synergia.matching.envelope_match(emitx,emity,current,gourmet)
    
    print " "
    print "  printed from  main"
    print " "
    print "sigma_x = %5.8f" % retval[0]
    print "sigma_xprime= %5.8f" % retval[1]
    print "r_x= %5.8f" % retval[2]
    print "sigma_y = %5.8f" % retval[3]
    print "sigma_yprime= %5.8f " % retval[4]
    print "r_y= %5.8f" % retval[5]
    
    

