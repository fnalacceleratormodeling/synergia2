#!/usr/bin/env python

from UberPkgpy import Apply_SpaceCharge_external as apply_space_charge_kick

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return retval
