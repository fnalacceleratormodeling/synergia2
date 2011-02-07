Overview
========

Synergia2 is a parallel particle-in-cell (PIC) beam dynamics framework.

Lattice
objects hold lists of unique instances of Lattice_elements, which are abstract
descriptions of the elements of a real accelerator. A Lattice_simulator object
contains the additional information necessary to perform simulations on Lattices.
Arbitrary string and floating-point data can be stored in the attributes of
Lattice_element objects.

A simulation consists of a series of Steps composed of the application of a
series of Operators to a macroparticle Bunch.
Operators are divided into independent and collective.
Only the former are required to be present. The application of
Independent_operators consists of one or more Independent_operations. Built-in
Independent_operations include chef_propagate and chef_map types. The type
of operation corresponding to a Lattice_element is determined by an
Operation_extractor_map, which uses information stored in the extractor_type
string attribute of the Lattice_element.

The simulation Steps are extracted by Stepper objects, of which several different
types exist corresponding to different stepping algorithms.

The simulation proceeds by passing the Stepper objects to a Propagator object,
which applies Steps to Bunch. The Propagator also performs and writes specified
Diagnostics on a per-step and per-turn basis.
