from .bunch import *
import numpy as np

def get_particles_numpy(self, group = ParticleGroup.regular):
    return np.array(self.get_host_particles(group), copy=False)

setattr(Bunch, 'get_particles_numpy', get_particles_numpy)


def calculate_mean(b):
    return np.array(Core_diagnostics.calculate_mean_ka(b))

def calculate_abs_mean(b):
    return np.array(Core_diagnostics.calculate_abs_mean_ka(b))

def calculate_std(b, mean):
    #mean = Core_diagnostics.calculate_mean_ka(b)
    return np.array(Core_diagnostics.calculate_std_ka(b, mean))

setattr(Core_diagnostics, 'calculate_mean', staticmethod(calculate_mean))
setattr(Core_diagnostics, 'calculate_abs_mean', staticmethod(calculate_abs_mean))
setattr(Core_diagnostics, 'calculate_std', staticmethod(calculate_std))
