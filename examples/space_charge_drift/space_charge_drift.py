#!/usr/bin/env python

# space charge drift test

import sys
import os

import synergia
import scipy
from scipy import interpolate
import numpy as np

from space_charge_drift_options import opts

mp = synergia.foundation.pconstants.mp
e = synergia.foundation.pconstants.e # [C]
c = synergia.foundation.pconstants.c # m/s
rp = synergia.foundation.pconstants.rp # [m]
mp = synergia.foundation.pconstants.mp # [GeV/c^2]
#  characteristic current I = 4 pi epsilon_0 m c^3/q
I0 = c * e/rp
print "I0: ", I0
#########################
#
#    Analytic calculations thanks to Nathan Cook and Chris Hall of RadiaSoft

##########################
# ## Analytical Comparison
##########################


# In[19]:

def calc_perveance(I,ref,cn=0):
    '''Calculate the perveance for a proton beam of a given current and particle energy.

    Arguments
        - I - current in A
        - ref - the reference particle for extracting beta and gamma

        - (optional) charge neutralization factor - default 0
    '''

    beta = ref.get_beta()
    gamma = ref.get_gamma()

    return (I/I0)*(2/beta**3)*(1/gamma**3)


# In[20]:

#Introduce numerical integrators

#2nd Order RK - Ralston Method
def Ralston(r,z,h,f):
    k1 = h*f(r)
    return 0.25*k1 + 0.75*h*f(r+(2/3)*k1)

#4th Order Runge-Kutta
def RungeKutta4(r,z,h,f):
    k1 = f(r)
    k2 = f(r + (h/2)*k1)
    k3 = f(r + (h/2)*k2)
    k4 = f(r + h*k3)
    return h/6*(k1 + 2*k2 +2*k3 + k4)

#function here, which is a function of r and z
def rprime(K,emit,r0,rp0,rm):
    '''

    Returns the slope of the beam envelope (dr/dz) for a given value of emittance,rm, K, and initial conditions.

    This equation follows from Reisier.

    Arguments:

        - r - beam radius (or RMS)
        - K - perveance
        - emit - geometric emittance
        - r0 - initial envelope radius (or RMS)
        - rp0 - initial slope of envelope (or RMS)

    '''

    first = rp0**2 #first term
    second = (emit**2)*((1./r0**2)-(1./rm**2)) #second term
    third = 2*K* np.log(rm/r0) / 4

    return np.sqrt(first + second + third)


# In[21]:cu

def calculate_expansion(current, reference_particle,r0,rp0,emit,N,zf):

    '''Evaluate the expansion of a KV beam envelope in a drift along z-axis, begining at z = 0.

    Arguments:
        - current - beam current in A
        - reference_particle - synergia object for bunch/lattice reference particle
        - r0 - initial envelope value (provide RMS for RMS expansion, a for envelope expansion, etc.)
        - rp0 - initial slope of envelope (must be non-zero, but calculation is not sensitive to small values)

        - (optional) emit - geometric emittance of beam - default 2.05721258396*1.e-6 (for 0.3 mm-mrad KV beam)
        - (optional) N - number of steps for integration - default 1000
        - (optional) zf - final z value (e.g. length of expansion) - default 50.0

    '''

    z0 = 0.0 #start
    ss = (zf-z0)/N #step size

    zpoints = np.linspace(0.0, zf, num=N) #define z values
    rpoints = [] #empty array for r values

    #calculate perveance
    Kp = calc_perveance(current, reference_particle)

    #x is r
    #z is t (what we step up)
    #f is our function describing the relationship between r and z
    f = lambda r: rprime(Kp,emit,r0,rp0,r)

    r,z,dz = r0,z0,ss
    points = []
    while z < zf:
        points.append((z,r))
        #z, r = z+dz, r + Ralston(r,z,dz,f) #incremement
        z, r = z+dz, r + RungeKutta4(r,z,dz,f) #incremement

    return scipy.array(points)




#########################
# create lattice of a single drift
lattice = synergia.lattice.Lattice("channel", synergia.lattice.Mad8_adaptor_map())
drift = synergia.lattice.Lattice_element("drift", "drift")
drift.set_double_attribute("l", opts.driftlength)
lattice.append(drift)

print "Lattice length: ", lattice.get_length()

#########################
# create the reference particle
reference_particle = synergia.foundation.Reference_particle(1, mp, opts.ke+mp)
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()
lattice.set_reference_particle(reference_particle)
print "Beam energy: ", reference_particle.get_total_energy()
print "Beam momentum: ", reference_particle.get_momentum()
print "Beam gamma: ", reference_particle.get_gamma()
print "Beam beta: ", reference_particle.get_beta()

#########################
# emittance from requested normalized emittance
emit = opts.nemit/(beta*gamma)
print "Beam normalized emittance: ", opts.nemit
print "Beam geometric emittance: ", emit
#########################
# calculations towards generating the bunch
#
# I = e * N/L * v
current = opts.current
real_particles = current*opts.blen/(e * beta * c)
#real_particles = opts.real_particles
#current = opts.real_particles * e * beta * c/opts.blen
print "Bunch length: ", opts.blen
print "Beam current: ", current
print "Beam bunch charge [e]: ", real_particles
commxx = synergia.utils.Commxx()
bunch = synergia.bunch.Bunch(reference_particle, opts.macroparticles, real_particles, commxx)
dist = synergia.foundation.Random_distribution(opts.seed, commxx)
# populate KV with same emittance x/y, uniform in cdt, monochromatic
synergia.bunch.populate_transverse_KV_GaussLong(dist, bunch, emit, 0.0, opts.betax, emit, 0.0, opts.betay, opts.blen/beta, 0.0)
local_particles = bunch.get_local_particles()
# center the distribution
means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
for i in range(5):
    local_particles[:, i] -= means[i]
diag_full2 = synergia.bunch.Diagnostics_full2("foo.h5")
diag_full2.set_bunch(bunch)
diag_full2.update()
means = diag_full2.get_mean()
stds = diag_full2.get_std()
mom2 = diag_full2.get_mom2()
emitx = diag_full2.get_emitx()
emity = diag_full2.get_emity()
del diag_full2
print "Bunch initial means: ", means
print "Bunch initial stds: ", stds
print "Bunch initial emitx: ", emitx
print "Bunch initial emity: ", emity

# dump bunch particles
diag_particles = synergia.bunch.Diagnostics_particles("kv_particles.h5")
diag_particles.set_bunch(bunch)
diag_particles.update_and_write()

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
bunch_simulator.add_per_step(synergia.bunch.Diagnostics_full2("full.h5"))


if opts.solver == "2d-openhockney":
    grid = [opts.gridx, opts.gridy, opts.gridz]
    coll_operator = synergia.collective.Space_charge_2d_open_hockney(commxx, grid)
    print "Space charge using 2d open hockney, grid: ", grid
elif opts.solver == "2d-kv":
    print "Space charge using 2d KV"
    coll_operator = synergia.collective.Space_charge_2d_kv()
else:
    raise RuntimeError, "Unknown solver: %s"%opts.solver

stepper = synergia.simulation.Split_operator_stepper(lattice, 1, coll_operator, opts.steps)

propagator = synergia.simulation.Propagator(stepper)
propagator.set_checkpoint_period(100)
propagator.propagate(bunch_simulator, opts.turns, 0, opts.verbosity)

final_means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
final_stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
print "Bunch final means: ", final_means
print "Bunch final stds: ", final_stds

del propagator
del stepper
del bunch_simulator
del bunch

# this is [ sqrt(E( (x_i + d*x'_i)**2)) - sqrt(E(x_i**2)) ]/d as the first guess for rprime 
r0 = stds[0]
rp0 = mom2[0,1]/scipy.sqrt(mom2[0,0])

calc_expansion = calculate_expansion(current, reference_particle, r0, rp0, emit, 1000000, opts.turns*lattice.get_length())

scipy.save("calc_expansion", calc_expansion)

fractional_diff = (final_stds[0] - calc_expansion[-1,1])/final_stds[0]
print "Fraction difference at end of channel (synergia - envelope)/synergia: ",  fractional_diff

if opts.plot:
    import matplotlib.pyplot as plt
    import tables

    h5 = tables.openFile("full.h5")
    stds = h5.root.std.read()
    s = h5.root.s.read()
    h5.close()

    plt.plot()
    plt.title("space charge expansion, current = %f [A]"%current)
    plt.plot(calc_expansion[:,0], calc_expansion[:,1], label='envelope equation')
    plt.plot(s, stds[0, :], label='Synergia solver=%s'%opts.solver)
    plt.legend(loc='best')
    plt.xlabel('s [m]')
    plt.ylabel("std x")
    plt.show()

