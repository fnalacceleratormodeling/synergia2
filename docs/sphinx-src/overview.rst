Overview
========

Synergia2 is a particle-in-cell (PIC) beam dynamics framework. Our aim is to simulate complex
problems in the realm of accelerator physics where analytical solutions, such as those obtained by
beam envelope equations, are not accurate enough.  Collective effects induced by space charge or
electron clouds, as well non-trivial aspects of single particle  optics, can be simulated
accurately and jointly in timely manner. Our ultimate goal is to support upgrades of existing high
energy accelerators and basic designs of new machines, where collective effects are important.

The Synergia2 framework is build on previous attempts and experiences in dealing with
challenging simulations of high intensity beams, such as those delivered by the Fermilab
booster.  As such, it became quickly an important integration task, where
algorithms used in single particle optics  have to be combined  with  collective effects.   As
such, Synergia2 is not delivered as a pre-build executable that runs on a desktop.  Synergia2 is
both a framework and tool-kit, supported by open-source code utilities, interfaces and
libraries target for specific problems.

Our aim is to predict the motion  of high energy particles in a beam (bunched or continuous), 
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

Familiarity with some basics of high energy accelerator_ science is strongly advised.  Using
Synergia2 does not require advanced knowledge of Object Oriented (OO) methods, nor
message-passing softwares.   Instead of filling a rigid set of input parameters, the Synergia2
user will compose (or better, modify) a  Python script of a simple C++ to describe the problem
at hands.  Such scripts can easily and seamlessly be upgraded to run on clusters.  Users are
invited to learn Synergia2 via the successful methodology of "learning by examples".  That is,
find the most relevant examples delivered in the package, run them, and then adapt.

As such, Synergia2 is also extensible:  Should a new collective effect occur, new modules and
be composed using the existing infra-structure.

We now review the basic concepts and elements of Synegia2. 

Lattice
-------
A Lattice is a computer model of a real accelerator. 
Objects (such a representaiton of a dipole magent, quadrupole, drift region, etc..) 
in this lattice hold lists of unique instances of ``Lattice_elements``. 
This list contains placement information for such devices.
A ``Lattice_simulator`` object contains the additional information necessary 
to perform simulations on Lattices.
Arbitrary string and floating-point data can be stored in the attributes of
Lattice_element objects.
Such a Lattice can be implemented by coding a Python or C++ program that creates 
and assembles the accelerator model. A more convenient way is to rely on our MAD_ input
file parser. More precisely, the MAD8 input files are fully supported [#]_.

A Synergia Run
--------------
A simulation consists of a series of Steps composed of the application of a
series of Operators to a macroparticle Bunch.
Operators are divided into independent and collective. By independent, we mean that each
macroparticle motion is treated independently, following the laws of single particle optics
in the accelerator.  Only the former are required to be present. Collective operators 
allows for a self-consistent simulation of space charge wake field effects.

A Synergia Bunch
----------------
A Bunch_ is an ensemble of macroparticles. As a typical bunch in a real accelerator 
can contain 1E11 to 1E13 particles, simulating the motion of all these particles 
is both prohibitively expensive and unnecessary. As E.M. fields do add linearly, 
a sub-ensemble ofparticles located in a small region of physical space will create and and be
sbjected to E.M. field that is equal to the one created (or seen)  
by one single particle with an electric charge commensurate to the sum of the
charge found in that region.  This single particles is refered as a macroparticle. 

Single Particle Optics
----------------------
The simulation of a single macroparticle in the accelerator is 
achieved via the use of ``Independent_operators``.
Such entities may consist of one or more changes in the position 
of the macroparticle in 6D phase space.  The tool used in Synergia for such
operations is CHEF_. 
Built-in Independent_operations include ``chef_propagate`` and ``chef_map types``. The type
of operation corresponding to a ``Lattice_element`` is determined by an
``Operation_extractor_map``, which uses information stored in the ``extractor_type``
string attribute of the ``Lattice_element``.


Synergia Steps
--------------
A Synergia simulation is performed in the time domain. That is, a Bunch_ propagates 
through the machine in discretes time steps. Such actions are managed by Stepper_ objects.
Informations about what occured in a give step can be accessed from such Stepper objects, 
of which several different types exist corresponding to different stepping algorithms.

The simulation proceeds by passing the Stepper objects to a Propagator object,
which applies Steps to Bunch_. The Propagator also performs and writes specified
Diagnostics on a per-step and per-turn basis.

.. _MAD: http://mad.home.cern.ch/mad/
.. _Bunch: ./bunch.html
.. _Stepper: ./Stepper.html
.. _CHEF: ./CHEF.html
.. _accelerator: http://uspas.fnal.gov/

.. rubric: Footnotes

.. [#]   We plan to support MADX files as well in the coming futur.


 
