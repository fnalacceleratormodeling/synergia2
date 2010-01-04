#!/usr/bin/env python

from workflow import Options,Job_manager
import workflow
import sys

if ( __name__ == '__main__'):
    myopts = Options("circular")
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",1,"map order",int)
    myopts.add("emittance95_over_pi",20,"95\% emittance in units of pi * mrad*mm",float)
    # longitudinal beta is 143.6
    myopts.add("dpop",3.482e-4,"(delta p)/p RMS width",float)
    #~ myopts.add("bunchlen", 0.05, "RMS bunchs length (z width) [m]", float)
    myopts.add("bunchlen", 40.0, "RMS bunch length (z width) [nanoseconds]", float)
    myopts.add("dpopoffset", 0.0, "offset in dpop", float)
    myopts.add("kicks",1,"kicks per line",int)
    myopts.add("turns",2,"number of turns",int)
    #~ myopts.add("latticefile","Debunch_modified.lat","",str)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",64,"",int)
    myopts.add("space_charge",1,"",int)
    myopts.add("impedance",0,"",int)
    myopts.add("energy",8.87710994,"total energy, default taken from Debunch_modified.lat",float)
    myopts.add("partpercell",1,"",float)
    #~ myopts.add("current",0.5,"current",float)
    myopts.add("realnum",1.2e13,"number of real particles per bunch",float)
    myopts.add("solver","3d","solver",str)
    myopts.add("aperture",0.05,"aperture radius in m",float)
    myopts.add("numtrack",1000,"number of particles to track",int)
    myopts.add("solver","3d","solver (3d or 2d)",str)
    myopts.add("rampturns",200,"sextupole ramping turns",int)
    myopts.add("periodic",0,"longitudinal periodic boundary conditions",int)
    myopts.add("transverse",0,"use longitudinally uniform beam",int)
    myopts.add("tuneh",9.63,"horizontal fractional tune",float)
    myopts.add("tunev",9.75,"vertical fractional tune",float)
    
    myopts.add_suboptions(workflow.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = Job_manager("mu2e.py",sys.argv,myopts,
                                      ["Debunch_modified.lat","dgourmet.py","debuncher.so"])
