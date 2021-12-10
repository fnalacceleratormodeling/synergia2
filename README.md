# Synergia

![](wiki/animation_synergia.gif)

`Synergia` is a accelerator modeling and simulation package developped at Fermilab.

- [Build Synergia](#markdown-header-build-instructions)
- [Writing a Synergia Simulation](#markdown-header-writing-a-synergia-simulation)
    - [A FODO Example in Python](#markdown-header-a-fodo-example-in-python)
    - [Explanation](#markdown-header-explanation)
        - [Lattice](#markdown-header-lattice)
        - [Bunch](#markdown-header-bunch)
        - [Populate Bunch Data](#markdown-header-populate-bunch-data)
        - [Bunch_simulator](#markdown-header-bunch_simulator)
        - [Diagnostics](#markdown-header-diagnostics)
        - [Propagator and Space Charge Operators](#markdown-header-prpagator-and-space-charge-operators)
        - [Propagate](#markdown-header-propagate)
        - [Run the Script](#markdown-header-run-the-script)
- [Advanced Topics](#markdown-header-advanced-topics)
- [API References](#markdown-header-api-and-class-references)


## Build Instructions 

Please see [this page](wiki/build.md) for build instructions

## Writing a Synergia Simulation

Here is a simple example to give an idea of how to compose a simulation with Synergia3.

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
    synergia.simulation.checkpoint_save(propagator, sim)

def main():

    try:
        run()
    except:
        raise RuntimeError("Failure to launch fodo.run")

main()
```

## Explanation

In the above example, a bunch of 1 million particles is created and populated with
the given mean and covariances. This bunch will propagate through a simple
focusing-drift-defocusing-drift, or `FODO` lattice for 10 turns. At each turn, a 
full particle diagnostic will be performed and saved in the `diag_full.h5`
file. After 10 turns, it writes a checkpoint save of the current state of the
simulation so it can resume later.

### Lattice

`Lattice` describes the structure of the accelerator complex that we are simulating. 
A `Lattice` object can be created by reading in a Mad8/MadX file or a string as shown 
in the above example. E.g.,

```python
# construct a MadX reader object
reader = synergia.lattice.MadX_reader()

# read in the madx from a string
reader.parse(madx_string)

# or can read in from a madx file
reader.parse_file(madx_file)

# finally, extract the lattice from the named sequence or line
lattice = reader.get_lattice(sequence_name)
```

Alternatively, a lattice object can also be created programmatically by appending
`Lattice_element` objects to the lattice manually,

```python
# lattice object needs a reference particle
ref_part = synergia.foundation.Reference_particle(charge, mass, energy)

# create a lattice object
lattice = synergia.lattice.Lattice(name, ref_part)

# create lattice elements
f = synergia.lattice.Lattice_element(type='quadrupole', name='f')
f.set_double_attribute('k1', 0.01)

o = synergia.lattice.Lattice_element(type='drift', name='o')
o.set_double_attribute('l', 1.0)

# and append elements to the lattice
lattice.append(f)
lattice.append(o)
```

### Bunch

`Bunch` is the container object for the particle data in Synergia. To create a `Bunch` 
object you need a `Reference_particle` for the bunch, the number of 
`macro_particles`, and the number of `real_particles`,

```python
bunch = synergia.bunch.Bunch(reference_particle, total_num, real_num, comm = Commxx())
```

where `total_num` is the number of macroparticles to be simulated in the bunch, 
`real_num` is the number of real particles of the bunch. `comm` is an optional argument
of an MPI communicator indicating the processes the bunch will span across. It is
defaulted to the `MPI_COMM_WORLD`.

Various methods are provided to access the particle data and properties of a bunch.
Commonly used ones are,

```python
# get the number of particles reside on the current rank
Bunch.get_local_num()

# get the total number of particles of this bunch
Bunch.get_total_num()

# retrieve the bunch reference particle
Bunch.get_reference_particle()

# retrieve the lattice design reference particle
Bunch.get_design_reference_particle()

# retrieve the particle data in a 2d numpy array [0:local_num, 0:6]
# in the second dimension, 0 - x, 1 - xp, 2 - y, 3 - yp, 4 - cdt, 5 - dpop, 6 - id
Bunch.get_host_particles_numpy()
```

When running a simulation on compute accelerators such as GPUs, it is crucial to know
that the particle data resides on both the device memory (such as GPU VRAM) and host
memory. All the propagation and computation are performed on the array of device 
memory. But user can only access the particle data on host memory. Therefore, a sync
operation is required if user seeks to read or modify the particle data during
propagation.

```python
# sync the particle data from host to device
Bunch.checkin_particles()

# sync the particle data from device to host
Bunch.checkout_particles()
```

### Populate Bunch Data

A newly created bunch has its all particle coordinates initialized to `0.0`. User may
choose to populate the particle coordinates by directly accessing the particle data
as shown in previous section. 

`Synergia` also provides a range of bunch populate methods to initialize the particle
to certain distributions. For example, the following one populate a bunch with given 
means and covariances,

```python
# mean and covariance matrix. The covariances matrix is 6x6
bunch_means = np.zero(6, dtype='d')
bunch_covariances = np.array([[...], [...], [...], [...], [...], [...]])

# use the provided random number generator
dist = synergia.foundation.PCG_random_distribution(seed)

# populate a Gaussian bunch
synergia.bunch.populate_6d(dist, bunch, bunch_meas, bunch_covariances)
```

A bunch can also be initialized from reading a generated particle file,

```python
bunch.read_file("turn_particles_0000.h5")
```

### Bunch_simulator

In order to propagate the generated bunch through a lattice, we need to do a little
extra work -- creating a `Bunch_simulator` object. A `Bunch_simulator` is the container
for holding one or multiple bunches that will be sent to the lattice for propagation.
It is also the object for registering diagnostics actions with the bunches.

A `Bunch_simulator` can be created from an existing bunch. It is however more 
convenient and optimal to create a `Bunch_simulator` with bunches in it. For the reason
that when running simulations with multiple bunches on multiple MPI processes the
built-in `construct()` method will do the decomposition and distribute bunches across
all processes in an optimal way.

```python
# reference particle, number of macro particles, and number of real particles
reference_particle = synergia.foundation.Reference_particle(charge, mass, energy)
macroparticles = 1024
realparticles = 1e13

# create the Bunch_simulator
simulator = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
    reference_particle, macroparticles, realparticles)

# the bunch object can be access from the simulator for further operations
bunch = simulator.get_bunch()
```

### Diagnostics

During the simulation, it is often critical to monitor the shape and properties of the
particle beam. In synergia it is achieved by registering various beam diagnostic
objects to the `Bunch_simulator`. `Synergia` provides common diagnostics such as
getting the statistics of a bunch, or writing out coordinates for selected particles, 
etc.

```python
# register a bunch statistic diagnostic performing at the end of every turn
diag_full2 = synergia.bunch.Diagnostics_full2("diag_full2.h5")
simulator.reg_diag_per_turn(diag_full2)

# write out the state of first 100 particles every other turn
diag_particles = synergia.bunch.Diagnostics_particles("diag_particles.h5", 100)
simulator.reg_diag_per_turn(diag_particles, period=2)
```

`Synergia` also allows user to register diagnostics on a per-step, or at a specific
lattice element by calling the following methods,

```python
simulator.reg_diag_per_step(diag, period)
simulator.reg_diag_at_element(diag, element)
```

### Propagator and Space Charge Operators

`Synergia` can simulate the beam dynamics with full three-dimensional space charge 
effect. It implements the split-operator technique for the space charge operators.
The object that employs this split-operator is the `Propagator`.

To create a `Propagator` we would need to prepare a machine lattice from the
`Lattice` object (as shown above), and a `Stepper` object to split the lattice and
insert the space charge operator at proper positions.

```python
# prepare the space charge operator
sc_ops = synergia.collective.Space_charge_3d_open_hockney_options(gridx, gridy, gridz)
sc_ops.comm_group_size = 1

# we will use the split operator stepper for the propagator
stepper = synergia.simulation.Split_operator_stepper(sc_ops, num_steps)

# finally we may create the propagator with lattice and stepper
propagator = synergia.simulation.Propagator(lattice, stepper)
```

### Propagate

With both the `Propagator` and `Bunch_simulator` objects ready, we may instruct the 
propagator to start the propagation,

```python
# a log object for writing out the system messages during simulation
logger = synergia.utils.parallel_utils.Logger(
    0, synergia.utils.parallel_utils.LoggerV.INFO_TURN)

# propagate for number of turns
propagator.propagate(simulator, logger, turns)
```

### Run the script

Now we have our first `Synergia` simulation script. Run the simulation with,

```
python fodo.py
```

It should produce outputs like,

```
Propagator: starting turn 1, final turn 10

Propagator: turn 1/inf., time = 0.278s, macroparticles = (1048576) / ()
Propagator: turn 2/inf., time = 0.247s, macroparticles = (1048576) / ()
Propagator: turn 3/inf., time = 0.300s, macroparticles = (1048576) / ()
Propagator: turn 4/inf., time = 0.351s, macroparticles = (1048576) / ()
Propagator: turn 5/inf., time = 0.308s, macroparticles = (1048576) / ()
Propagator: turn 6/inf., time = 0.275s, macroparticles = (1048576) / ()
Propagator: turn 7/inf., time = 0.272s, macroparticles = (1048576) / ()
Propagator: turn 8/inf., time = 0.259s, macroparticles = (1048576) / ()
Propagator: turn 9/inf., time = 0.249s, macroparticles = (1048576) / ()
Propagator: turn 10/inf., time = 0.263s, macroparticles = (1048576) / ()
Propagator: maximum number of turns reached
Propagator: total time = 2.940s
```


## Advanced Topics

Here some advanced topics for Synergia simulation.

### Lattice_simulator

TBA

### Spectator Particles

TBA

### Propagate Actions

TBA

### Diagnostics and Custom Diagnostics

TBA

### Longitudinal Boundaries

TBA

### Bunch Injection

TBA

### Bunch Train and Multiple Bunch Simulation

TBA

## API and Class References

TBA






