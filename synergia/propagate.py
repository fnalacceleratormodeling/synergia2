#!/usr/bin/env python

import s2_fish
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
def propagate(s0,gourmet,bunch_in,diagnostics_in,grid_dim,quiet=1,
    use_s2_fish=False, use_impact=False, use_none=False,
    use_s2_fish_cylindrical=False,
    pgrid=None,field=None,cgrid=None,use_gauss=False,
    periodic=False, aperture=None, radius=None,
    space_charge=True,impedance=False,
    pipe_radiusx=None,pipe_radiusy=None,
    pipe_conduct=None,bunch_spacing=None,
    tracker=None,track_period_steps=None,
              transverse=False):

    bunches = listify(bunch_in)
    diagnosticss = listify(diagnostics_in)
    
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
                #~ pardebug("start space charge\n")
                if use_s2_fish:
                    s2_fish.apply_space_charge_kick(grid_dim,None,None, bunches, 2*tau,
                        periodic=periodic,aperture=aperture,space_charge=space_charge,
                        impedance=impedance,pipe_radiusx=pipe_radiusx,
                        pipe_radiusy=pipe_radiusy,pipe_conduct=pipe_conduct,
                        bunch_spacing=bunch_spacing,transverse=transverse)
                elif use_s2_fish_cylindrical:
                    s2_fish.apply_cylindrical_space_charge_kick(grid_dim,
                        radius,bunch,2*tau,aperture=aperture,space_charge=space_charge,
                        impedance=impedance,impedance_pipe_radiusx=pipe_radiusx,
                        impedance_pipe_radiusy=pipe_radiusy,
                        pipe_conduct=pipe_conduct,
                        bunch_spacing=bunch_spacing)
                elif use_impact:
                    if not have_impact:
                        raise RuntimeError, \
                            "propagate with use_impact=True requires a working impact module"
                    if impedance:
                        raise RuntimeError,\
                            "impact solvers cannot currently be combined with impedance"
                    if ((pgrid == None) or (field == None) or (cgrid == None)):
                        raise RuntimeError, \
                            "propagate with use_impact=True requires pgrid, field and cgrid to be specified"
                    impact.apply_space_charge_kick(
                        bunch.get_beambunch(),
                        pgrid.get_pgrid2d(),
                        field.get_fieldquant(),
                        field.get_compdom(),
                        field.get_period_length(),
                        cgrid.get_bc_num(),
                        field.get_pipe_radius(),
                        tau, 0, bunch.get_scaling_frequency(),0)
                elif use_gauss:
                    s2_fish.apply_BasErs_space_charge_kick(bunch, 2*tau,
                        space_charge=space_charge,
                        impedance=impedance,impedance_pipe_radiusx=pipe_radiusx,
                        impedance_pipe_radiusy=pipe_radiusy,
                        pipe_conduct=pipe_conduct,
                        bunch_spacing=bunch_spacing)
                elif use_none:
                    pass
                else:
                    raise RuntimeError, \
                        "propagate requires one of use_s2_fish, use_impact or use_none to be True"
                #~ pardebug("end space charge\n")
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
