Overview
========

Synergia2 is a particle-in-cell (PIC) beam dynamics framework. Our aim is to simulate complex
problems in the realm of accelerator physics, where analytical solutions, such as those
obtained by beam envelope equations, are not accurate enough.  Collective effects, such as
those induced by space charge or electron clouds, as well non-trivial aspects of single
particle  optics, can be simulated accurately in timely manner. Our ultimate goal is to support
upgrades of existing high energy accelerators and basic designs of new machines, where
collective effects are important.

The Synergia2 framework is build on previous attempts and experiences in dealing with
challenging simulations of high intensity beams, such as those delivered by the Fermilab
booster.  As such, it became quickly a framework and an important integration task, where
algorithm used in single particle optics  have to be combined  with  collective effects.   As
such Synergia2 is not delivered as a pre-build executable that runs on a desktop.  Synergia2 is
both a framework and tool-kit, supported by open-source code utilities, interfaces and
libraries target for specific problems.

Our aims is to predict the motion  of high energy particles in a beam (bunched or continuous), 
in 6D phase space.  For instance, Fields are expressed on a 3D rectangular grid, and, at any
given time, both longitudinal and transverse motions are treated consistently.   However, 
Synergia2 can also address simpler problems, where the longitudinal dynamics is clearly
disconnected from the transverse one, without paying an overwhelming price for the 3D aspect.

Finally, an essential part of the Synergia2 framework is its support for parallelism.  The
Synergia2 architecture has been designed from the onset to run on medium size clusters (few
hundreds to thousands  core or processors) equipped with adequate networking capabilities.  (
Our ultimate aims is to run on the next generation of leadership class machine, despite the
tremendous difficulties in the optimization of tightly coupled, 6D, systems).  However,
libraries and procedures to build and run Synergia2 on single desktops are also fully
supported.

Familiarity with some basics of  high energy accelerator science is strongly advised.  Using
Synergia2 does not require advanced knowledge of Object Oriented (OO) methods, nor
message-passing softwares.   Instead of filling a rigid set of input parameters, the Synergia2
user will compose (or better, modify) a  Python script of a simple C++ to describe the problem
at hands.  Such scripts can easily and seamlessly be upgraded to run on clusters.  Users are
invited to learn Synergia2 via the successful methodology of ``learning by examples".  That is,
find the most relevant examples delivered in the package, run them, and then adapt.

As such, Synergia2 is also extensible:  Should a new collective effect occur, new modules and
be composed using the existing infra-structure.

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
