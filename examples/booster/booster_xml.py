#!/usr/bin/env python
import synergia
from booster_options import opts

#quad_correctors_H=[] 
#quad_correctors_V=[] 
#sextupole_correctors_H=[] 
#sextupole_correctors_V=[] 
lattice=synergia.lattice.Mad8_reader().get_lattice("machine", "booster_2012.lat")
lattice_length=lattice.get_length() 
reference_particle = lattice.get_reference_particle()
beta = lattice.get_reference_particle().get_beta()
harmon=opts.num_buckets
freq=harmon*beta*synergia.foundation.pconstants.c/lattice_length
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

#    names= elem.get_name()[0:4] 
#    if ((names=="ssxs")  and (elem.get_type()=="sextupole")):        
#      sextupole_correctors_H.append(elem) 
#      #print "name=",elem.get_name() 
#    if ((names=="ssxl")   and (elem.get_type()=="sextupole")):        
#      sextupole_correctors_V.append(elem) 
#    #print "name=",elem.get_name()
#    nameq= elem.get_name()[0:3] 
#    nameql= elem.get_name()[0:6]
#    if ((nameq=="qql")  and (elem.get_type()=="quadrupole") and (nameql != "qqlerr") ):                    
#      quad_correctors_H.append(elem)     

#    if ((nameq=="qqs")  and (elem.get_type()=="quadrupole") and (nameql != "qqserr")):                    
#      quad_correctors_V.append(elem)       


synergia.lattice.xml_save_lattice(lattice, "booster_lattice.xml")

