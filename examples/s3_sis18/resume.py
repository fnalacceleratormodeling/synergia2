
import synergia as syn

def resume():

    simlog = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.INFO_STEP)

    screen = syn.utils.parallel_utils.Logger(0, 
            syn.utils.parallel_utils.LoggerV.DEBUG)

    [prop, sim] = syn.simulation.checkpoint_load()
    prop.propagate(sim, simlog, 1)

    syn.utils.parallel_utils.simple_timer_print(screen)

print("resuming sis_18")
resume()

