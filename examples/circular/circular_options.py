#!/usr/bin/env python

import workflow

myopts = workflow.Options("circular")
myopts.add("transverse", 0, "longitudinally uniform beam", int)
myopts.add("maporder", 1, "map order", int)
myopts.add("emittance", 5.89533703303356e-07, "emittance", float)
# longitudinal beta is 143.6
myopts.add("dpop", 3.482e-4, "(delta p)/p RMS width", float)
myopts.add("bunchlen", 0.05, "RMS bunchs length (z width) [m]", float)
myopts.add("dpopoffset", 0.0, "offset in dpop (am: ! -duop)", float)
myopts.add("kicks", 32, "kicks per line", int)
myopts.add("turns", 10, "number of turns", int)
myopts.add("latticefile", "foborodobo_s.lat", "", str)
myopts.add("tgridnum", 16, "transverse grid cells", int)
myopts.add("lgridnum", 32, "", int)
myopts.add("xoffset", 0.0004, "transverse offset in x", float)
myopts.add("yoffset", -0.0002, "transverse offset in y", float)
myopts.add("xpoffset", 0, "offset in x-prime", float)
myopts.add("ypoffset", 0, "offset in y-prime", float)
#    myopts.add("zoffset",0,"offset in z", float)
#    myopts.add("xoffset",4.26e-4,"transverse offset in x",float)
#    myopts.add("yoffset",1.86e-4,"transverse offset in y",float)
myopts.add("zoffset", 0.1, "offset in z", float)
# myopts.add("zoffset",0.,"offset in z", float)
myopts.add("space_charge", 0, "", int)
myopts.add("impedance", 0, "", int)
myopts.add("pipe_symmetry", "circular", "", str) 
#  myopts.add("pipe_symmetry","x_parallel_plates","",str) 
myopts.add("energy", 100.004401675138, "", float)
myopts.add("partpercell", 1, "", float)
myopts.add("bunches", 1, "", int)
myopts.add("bunchnp", 1.0e11, "number of particles per bunch", float)
myopts.add("kick", "full", "kick type", str)
myopts.add("prev_turns", 1, "number of prev turns contributing to impedance", int)

job_mgr = workflow.Job_manager("circular.py", myopts,
                               [myopts.get("latticefile")])
