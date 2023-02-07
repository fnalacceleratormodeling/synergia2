
import numpy as np
import synergia as syn

def run():

    # logger
    screen = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.DEBUG)

    simlog = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.INFO_TURN)

    # lattice
    with open('booster_init_lattice.json', 'r') as f:
        jslattice = f.read()
    lattice = syn.lattice.Lattice.load_from_json(jslattice)
    print('lattice length: ', lattice.get_length())
    print('lattice energy: ', lattice.get_reference_particle().get_total_energy())
    syn.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    # calculate one turn map at given order (3 here)
    mapping = syn.simulation.Lattice_simulator.get_one_turn_map_o3(lattice)

    # save
    mapping.save_json("mapping.json")

    # load
    m2 = syn.foundation.TMapping_o3.load_json("mapping.json")

    # create the normal form object from one turn mapping directly
    ref = lattice.get_reference_particle()
    e0 = ref.get_total_energy()
    pc0 = ref.get_momentum()
    mass = ref.get_mass()

    nf = syn.foundation.NormalForm_o3(m2, e0, pc0, mass)

    # save nf
    nf.save_json("nf.json")

    # load
    nf2 = syn.foundation.NormalForm_o3.load_json("nf.json")

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




