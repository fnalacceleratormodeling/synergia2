#!/usr/bin/env python

import s2_fish
import impact
import chef_propagate

last_step_length = 0
def propagate(s0,gourmet,bunch,diagnostics,grid_dim,quiet=1,
    use_s2_fish=False, use_impact=False, use_none=False,
    pgrid=None,field=None,cgrid=None,use_gauss=False):
    s = s0
    global last_step_length
    first_action = 1
    for action in gourmet.get_actions():
        if action.is_mapping():
            action.get_data().apply(bunch.get_local_particles(),
                               bunch.get_num_particles_local())
            last_step_length = action.get_length()
            s += last_step_length
        elif action.is_synergia_action():
            if action.get_synergia_action() == "space charge endpoint":
                if not first_action:
                    diagnostics.add(s,bunch)
                    if not quiet:
                        print "finished space charge kick"
            elif action.get_synergia_action() == "space charge kick":
                tau = last_step_length
                if use_s2_fish:
                    s2_fish.apply_space_charge_kick(grid_dim,None,None, bunch, 2*tau)
                elif use_impact:
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
                    s2_fish.apply_BasErs_space_charge_kick(bunch, 2*tau)
                elif use_none:
                    pass
                else:
                    raise RuntimeError, \
                        "propagate requires one of use_s2_fish, use_impact or use_none to be True"
            elif action.get_synergia_action() == "rfcavity1" or \
                action.get_synergia_action() == "rfcavity2":
                element = action.get_data()
                u_in = gourmet.get_u(action.get_initial_energy())
                u_out = gourmet.get_u(action.get_final_energy())
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
