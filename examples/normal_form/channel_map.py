
import numpy as np
import synergia as syn

def run():

    # logger
    screen = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.DEBUG)

    simlog = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.INFO_TURN)

    # lattice
    reader = syn.lattice.MadX_reader()
    lattice = reader.get_lattice("fodo", "channel.madx")
    syn.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    # calculate one turn map at given order (3 here)
    mapping = syn.simulation.Lattice_simulator.get_one_turn_map_o3(lattice)

    # convert the mapping to a json object
    mapping_json = mapping.to_json()
    print(mapping_json)

    # or iterate through the components and fields
    for comp in range(6):
        trigon = mapping.component(comp)
        print("\n\nComponent {}".format(comp))

        for pwr in range(trigon.power()+1):
           for idx in range(trigon.count(pwr)): 
               idx_to_exp = "syn.foundation.Trigon_index_to_exp_o{}(idx)".format(pwr)

               print('power = {}, exp = {}, term = {}'.format( 
                    pwr, eval(idx_to_exp), trigon.get_term(pwr, idx)))



print("channel_map")
run()




