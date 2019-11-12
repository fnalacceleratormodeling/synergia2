from bunch import *
import numpy as np

def get_particles_numpy(self):
    return np.array(self.get_host_particles(), copy=False)

setattr(Bunch, 'get_particles_numpy', get_particles_numpy)
