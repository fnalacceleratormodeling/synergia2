#!/usr/bin/env python
import sys
import s2_fish
from s2_fish.space_charge import *
have_impact = False
try:
    import impact
    have_impact = True
except ImportError:
    pass
import chef_propagate
import time

def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]

last_step_length = 0
def propagate(s0,gourmet,bunch_in,diagnostics_in, space_charge=None,impedance=None,
  aperture=None, bunch_spacing=None, track_period_steps=None,tracker=None, quiet=1):
  
    bunches = listify(bunch_in)
    diagnosticss = listify(diagnostics_in)
    #diagnosticss2 = listify(diagnostics2_in)
    
    if tracker and not track_period_steps:
        raise RuntimeError,\
            "propagate: when specifying a tracker, one must also specify track_period_steps"
    trackers=listify(tracker)

    if len(bunches) != len(diagnosticss):
        raise RuntimeError,\
            "propagate: len(bunches) must = len(diagnosticss)"
    s = s0
    steps = 0
    global last_step_length
    first_action = 1
    
    apply_fish_kick=False
    apply_impact_kick=False    
    if (space_charge) or (impedance):
        apply_fish_kick=True
        if space_charge:
            if (space_charge.get_solver()=="impact"):
                if impedance:
                       raise RuntimeError,\
                             "impact solvers cannot currently be combined with impedance"
                apply_impact_kick=True
                apply_fish_kick=False
    
     
   
    for action in gourmet.get_actions():
        if action.is_mapping():
            for bunch in bunches:
                action.get_data().apply(bunch.get_local_particles(),
                                   bunch.get_num_particles_local())
            last_step_length = action.get_length()
            s += last_step_length
        elif action.is_synergia_action():
            if action.get_synergia_action() == "space charge endpoint":
                if aperture:
                    for mbunch in bunches:
                        mbs = mbunch.get_store()
                        s2_fish.constraints.apply_circular_aperture(mbs,aperture)
                        mbunch.local_num = mbs.local_num
                        mbunch.total_num = mbs.total_num
                if not first_action:
                    if tracker:
                         if steps % track_period_steps == 0:
                             for (bunch,tracker) in zip(bunches,trackers):
                                 tracker.add(bunch,s)			
                    for (diagnostics,bunch) in zip(diagnosticss,bunches):
                            diagnostics.add(s,bunch)
                    if not quiet:
                        print "finished space charge kick"
                steps += 1
            elif action.get_synergia_action() == "space charge kick":
                tau = last_step_length
                if apply_fish_kick: 
		#~ 	s2_fish.apply_kick is in fish_kick.py
                    s2_fish.apply_kick(bunches, 2.*tau,space_charge=space_charge,
                        impedance=impedance,aperture=aperture) 
                elif apply_impact_kick:
                    apply_space_charge_kick(bunch,space_charge,tau)	#it skips aperture and periodic on the bunch               
		#for (diagnostics,bunch) in zip(diagnosticss2,bunches):
                        #diagnostics.add(s,bunch)	
            elif action.get_synergia_action() == "rfcavity1" or \
                action.get_synergia_action() == "rfcavity2":
                #~ pardebug("rfcavity\n")
                element = action.get_data()
                u_in = gourmet.get_u(action.get_initial_energy())
                u_out = gourmet.get_u(action.get_final_energy())
                for bunch in bunches:
                    chef_propagate.chef_propagate(
                        bunch.get_local_particles(), bunch.get_num_particles_local(),
                        element, action.get_initial_energy(), gourmet.particle,
                        u_in, u_out)
            else:
                print "unknown action: '%s'" % \
                      action.get_synergia_action()
        else:
            print "action",action.get_type(),"unknown"
        first_action = 0
   
    return s
