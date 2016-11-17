#!/usr/bin/env python
# -*- coding: utf-8 -*-
import synergia
from booster_options import opts
import numpy as np

def convert_orbum_to_drifts(orig_lattice):
    adaptor_map=synergia.lattice.MadX_adaptor_map()
    lattice = synergia.lattice.Lattice("rrnova",adaptor_map)
    for elem in orig_lattice.get_elements():
        if elem.get_name()[:9] == "orbumpa_2":
            new_elem = synergia.lattice.Lattice_element("drift", elem.get_name())
            s_attributes = elem.get_string_attributes()
            d_attributes = elem.get_double_attributes()
            for s in s_attributes.keys():
                new_elem.set_string_attribute(s, s_attributes[s])
            for d in d_attributes.keys():
                new_elem.set_double_attribute(d, d_attributes[d])
            lattice.append(new_elem)
        else:
            lattice.append(elem)
    lattice.set_reference_particle(orig_lattice.get_reference_particle())
    return lattice

def convert_dogs_to_drifts(orig_lattice):
    adaptor_map=synergia.lattice.MadX_adaptor_map()
    lattice = synergia.lattice.Lattice("rrnova",adaptor_map)
    for elem in orig_lattice.get_elements():
        if elem.get_name()[:3] == "dog":
            new_elem = synergia.lattice.Lattice_element("drift", elem.get_name())
            s_attributes = elem.get_string_attributes()
            d_attributes = elem.get_double_attributes()
            for s in s_attributes.keys():
                new_elem.set_string_attribute(s, s_attributes[s])              
            for d in d_attributes.keys():                
                if (d=='l'):
                 new_elem.set_double_attribute(d, d_attributes[d])
            lattice.append(new_elem)
        else:
            lattice.append(elem)
    lattice.set_reference_particle(orig_lattice.get_reference_particle())
    return lattice


def convert_rbends_to_sbends(orig_lattice):
    adaptor_map=synergia.lattice.MadX_adaptor_map()
    lattice = synergia.lattice.Lattice("rrnova",adaptor_map)
    for elem in orig_lattice.get_elements():
        if elem.get_type() == "rbend":            
            new_elem = synergia.lattice.Lattice_element("sbend", elem.get_name())
            s_attributes = elem.get_string_attributes()
            d_attributes = elem.get_double_attributes()
            for s in s_attributes.keys():
                new_elem.set_string_attribute(s, s_attributes[s])
            for d in d_attributes.keys():
                new_elem.set_double_attribute(d, d_attributes[d])
            
            if ( (not elem.has_double_attribute("angle") )  or  \
             ( elem.has_double_attribute("angle") and (abs(elem.get_double_attribute("angle"))<5.e-9) )  ):
                   arclength = elem.get_double_attribute("l")
                   ang=0.
            else:
                 ang = elem.get_double_attribute("angle") 
                 length = elem.get_double_attribute("l")
                 arclength = ang*length/(2.0*np.sin(ang/2.0)) 
          
            new_elem.set_double_attribute("angle", ang)
            new_elem.set_double_attribute("l", arclength)
            new_elem.set_double_attribute("e1", ang/2.0)
            new_elem.set_double_attribute("e2", ang/2.0)
            lattice.append(new_elem)
            #print "*******************"
            #new_elem.print_()
            #print "old"
            #elem.print_()
            #print "*******************"
        else:
            lattice.append(elem)

    lattice.set_reference_particle(orig_lattice.get_reference_particle())
    return lattice

def convert_zero_angle_sbends_to_drifts(orig_lattice):
    adaptor_map=synergia.lattice.MadX_adaptor_map()
    lattice = synergia.lattice.Lattice("rrnova",adaptor_map)
    for elem in orig_lattice.get_elements():
        if elem.get_type() == "sbend":
                if  ((elem.has_double_attribute("angle")) and (abs(elem.get_double_attribute("angle"))>1e-8)):
                        lattice.append(elem)
                else:
                      new_elem = synergia.lattice.Lattice_element("drift", elem.get_name())
                      s_attributes = elem.get_string_attributes()
                      d_attributes = elem.get_double_attributes()
                      for s in s_attributes.keys():
                          new_elem.set_string_attribute(s, s_attributes[s])
                      for d in d_attributes.keys():
                        new_elem.set_double_attribute(d, d_attributes[d])
                      lattice.append(new_elem)
        else: 
          lattice.append(elem)   
    lattice.set_reference_particle(orig_lattice.get_reference_particle())
    return lattice

def convert_zero_angle_rbends_to_drifts(orig_lattice):
    adaptor_map=synergia.lattice.MadX_adaptor_map()
    lattice = synergia.lattice.Lattice("rrnova",adaptor_map)
    for elem in orig_lattice.get_elements():
        if elem.get_type() == "rbend":
                if  ((elem.has_double_attribute("angle")) and (abs(elem.get_double_attribute("angle"))>1e-8)):
                        lattice.append(elem)
                else:
                      new_elem = synergia.lattice.Lattice_element("drift", elem.get_name())
                      s_attributes = elem.get_string_attributes()
                      d_attributes = elem.get_double_attributes()
                      for s in s_attributes.keys():
                          new_elem.set_string_attribute(s, s_attributes[s])
                      for d in d_attributes.keys():
                        new_elem.set_double_attribute(d, d_attributes[d])
                      lattice.append(new_elem)
        else: 
          lattice.append(elem)   

    lattice.set_reference_particle(orig_lattice.get_reference_particle())
    return lattice



#quad_correctors_H=[] 
#quad_correctors_V=[] 
#sextupole_correctors_H=[] 
#sextupole_correctors_V=[] 
lattice = synergia.lattice.MadX_reader().get_lattice("booster", "booster_synergia.madx")
#lattice=convert_orbum_to_drifts(lattice)
#lattice=convert_zero_angle_sbends_to_drifts(lattice)
lattice=convert_zero_angle_rbends_to_drifts(lattice)
lattice=convert_rbends_to_sbends(lattice)
#lattice=convert_dogs_to_drifts(lattice)
lattice_length=lattice.get_length() 
reference_particle = lattice.get_reference_particle()
beta = lattice.get_reference_particle().get_beta()
harmon=opts.num_buckets
freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
print "latice length=",lattice_length
print "beta=",beta
print "gamma=",lattice.get_reference_particle().get_gamma()
print "initial frequency=",freq      
for elem in lattice.get_elements():
    if elem.get_name()[:4] == "sept":       
        elem.print_()   
    if elem.get_name()[:3] == "dog":
       #elem.set_double_attribute("e1", 0)
       #elem.set_double_attribute("e2", 0)
       elem.print_()  
    #if elem.get_name()[:5] == "muldg":
    #   elem.print_()

    if opts.chef_propagate:
      elem.set_string_attribute("extractor_type", "chef_propagate")
    elif  opts.chef_map:   
       elem.set_string_attribute("extractor_type", "chef_map")

    if elem.get_name()[:2] == "rf":
       elem.set_double_attribute("volt", opts.rf_voltage)
      # elem.set_double_attribute("freq",  freq)
       elem.set_double_attribute("lag",  0.) 
       print "element=",elem.print_()
   
    if opts.aperture:   
      name= elem.get_name()[0:4]
      if (name=="fmag"):
         #print "name=",elem.get_name()
        elem.set_string_attribute("aperture_type","rectangular")
        elem.set_double_attribute("rectangular_aperture_width", 0.2)
        elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureF"))
      elif (name=="dmag"):
        ## print "name=",elem.get_name()
        elem.set_string_attribute("aperture_type","rectangular")
        elem.set_double_attribute("rectangular_aperture_width", 0.2)
        elem.set_double_attribute("rectangular_aperture_height", 2*opts.get("apertureD"))
      else:
        elem.set_string_attribute("aperture_type","circular")
        elem.set_double_attribute("circular_aperture_radius", opts.get("apertureL"))   


    #print "element=",elem.print_()
##    names= elem.get_name()[0:4] 
##    if ((names=="ssxs")  and (elem.get_type()=="sextupole")):        
##      sextupole_correctors_H.append(elem) 
##      #print "name=",elem.get_name() 
##    if ((names=="ssxl")   and (elem.get_type()=="sextupole")):        
##      sextupole_correctors_V.append(elem) 
##    #print "name=",elem.get_name()
##    nameq= elem.get_name()[0:3] 
##    nameql= elem.get_name()[0:6]
##    if ((nameq=="qql")  and (elem.get_type()=="quadrupole") and (nameql != "qqlerr") ):                    
##      quad_correctors_H.append(elem)     

##    if ((nameq=="qqs")  and (elem.get_type()=="quadrupole") and (nameql != "qqserr")):                    
##      quad_correctors_V.append(elem)       


synergia.lattice.xml_save_lattice(lattice, "booster_test.xml")

