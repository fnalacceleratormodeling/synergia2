#!/usr/bin/env python

import s2_fish
import chef_propagate

last_step_length = 0
def propagate(s0,gourmet,bunch,diagnostics,grid_dim,quiet=1):
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
                diagnostics.add(s,bunch)
                if not first_action and not quiet:
                    print "finished space charge kick"
            elif action.get_synergia_action() == "space charge kick":
                tau = last_step_length
                s2_fish.apply_space_charge_kick(grid_dim,None,None, bunch, 2*tau)
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