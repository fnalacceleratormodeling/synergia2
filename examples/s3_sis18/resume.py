
import synergia as syn

def resume():

    simlog = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.INFO_STEP)

    [prop, sim] = syn.simulation.checkpoint_load()
    prop.propagate(sim, simlog, 1)

print("resuming sis_18")
resume()

