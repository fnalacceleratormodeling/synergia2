# Synergia3

Synergia is a accelerator modeling and simulation package developped at Fermilab.

- [Build Synergia](#build-instructions)
- [Basic Concepts](#basic-concepts)
  - [Bunch](#bunch)
  - [Lattice](#lattice)
  - [Bunch_simulator and Propagator](#bunch_simulator-and-propagator)
- [Examples](#examples-of-synergia-simulations)
  - [FODO in Python](#a-fodo-example-in-python)
  - [FODO in C++](#fodo-in-c++)


## Build Instructions

Please see [this page](wiki/build.md) for build instructions

## Basic Concepts

### Bunch

`Bunch` is the object which holds the particle data in Synergia. To create a `Bunch` object you need a `Reference_particle` for the bunch, the number of `macro_particles`, and the number of `real_particles`.

### Lattice

`Lattice` describes the structure of a acclerator complex, along with a lattice design reference particle. A `Lattice` object can be created by reading in a Mad8/MadX file. Or alternatively, it can also be constructed by adding `Lattice_element` objects through the `Lattice::append()` method.

### Bunch_simulator and Propagator

In order to simulate the dynamics of a bunch (or a series of bunches) of particles propagating through an accelerator lattice, you will need two additional objects, the `Bunch_simulator`, and the `Propagator`.

`Bunch_simulator` holds the bunch that is going to be propagated, along with all the diagnostics and propagate actions that happen during the propagation.

`Propagator` contains the lattice structure for the simulation, and optional operators (such as the space charge solver, or apertures) that are inserted in between the lattice elements.

Finally the propagation happens when the `Propagator::propagate()` method is called, with the `Bunch_simulator` object as one of the arguments - it propagate the bunch along the lattice that is stored in the propagator object.

## Examples of Synergia Simulations

Here are some examples to give an idea of how to compose a simulation with Synergia3.

### A FODO example in Python

```python
#!/usr/bin/env python3

import numpy as np
import synergia

macroparticles = 1048576
real_particles=2.94e12
turns = 10
gridx = 32
gridy = 32
gridz = 128

def get_lattice():

    fodo_madx = """
      beam, particle=proton,pc=3.0;

      o: drift, l=8.0;
      f: quadrupole, l=2.0, k1=0.071428571428571425;
      d: quadrupole, l=2.0, k1=-0.071428571428571425;

      fodo: sequence, l=20.0, refer=entry;
      fodo_1: f, at=0.0;
      fodo_2: o, at=2.0;
      fodo_3: d, at=10.0;
      fodo_4: o, at=12.0;
      endsequence;
    """

    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice('fodo')
    return lattice

def create_simulator(ref_part):

    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
            ref_part, macroparticles, real_particles)

    bunch = sim.get_bunch()

    bunch_means = np.zeros(6, dtype='d')
    bunch_covariances = np.array(
        [[3.0509743977035345e-05, 2.2014134466660509e-06, 0, 0, 0, 0],
         [2.2014134466660509e-06, 1.9161816525115869e-07, 0, 0, 0, 0],
         [0, 0, 7.5506914064526925e-06, -6.6846812465678249e-07, 0, 0],
         [0, 0, -6.6846812465678249e-07, 1.9161816525115867e-07, 0, 0],
         [0, 0, 0, 0, 0.00016427607645871527, 0],
         [0, 0, 0, 0, 0, 1e-08]])
    
    dist = synergia.foundation.PCG_random_distribution(1234567)

    synergia.bunch.populate_6d(dist, 
        bunch, 
        bunch_means,
        bunch_covariances)

    return sim

def create_propagator(lattice):
    sc_ops = synergia.collective.Space_charge_3d_open_hockney_options(gridx, gridy, gridz)
    sc_ops.comm_group_size = 1

    stepper = synergia.simulation.Split_operator_stepper_elements(sc_ops, 1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

def run():

    # logger
    screen = synergia.utils.parallel_utils.Logger(0, 
      synergia.utils.parallel_utils.LoggerV.INFO)

    simlog = synergia.utils.parallel_utils.Logger(0, 
      synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    # components
    lattice = get_lattice()
    sim = create_simulator(lattice.get_reference_particle())
    propagator = create_propagator(lattice)

    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag_full.h5")
    sim.reg_diag_per_turn(diag_full2)

    # propagate
    propagator.propagate(sim, simlog, turns)

    # save
    synergia.simulation.checkpoint_save(propagator, sim);

def main():

    try:
        run()
    except:
        raise RuntimeError("Failure to launch fodo.run")

main()
```

In the above example, a bunch of 1 million particles is created and populated with the given mean and covariances. This bunch will propagate through a simple focusing-drift-defocusing-drfit, or `FODO` lattice for a given number of turns. At each turn, a full particle diagnostics will be performed and saved in the `diag_full.h5` file. After 10 turns, it writes a checkpoint save of the current state of the simulation so it can be resumed later.


### FODO in C++

This is the same `FODO` example but written in C++.

