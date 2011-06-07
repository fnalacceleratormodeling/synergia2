Examples
========

fodo
----

This is the simplest example the Synergia authors could think off. Such a system has been described on various
occasions (see for instance the notes from the USPAS accelerator school course by Mike Syphers_ et al). Both
single particle optic and collective effect are demonstrated. This example runs in a few minutes on a single CPU
(or core) system. 

The lattice is of the simple Focus-Defocus (FoDo) type: one focusing quadrupole followed by a short drift, 
then one defocusing quadrupole of the same strength and ending with a drift of the same length. 


To get going, first, execute the setup shell script found in the top level directory where the Synergia package
was installed::

  % source .../Synergia/synergia2-refactor/setup.sh 
 

A. Via The Python Interface. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The lattice is defined in just a few lines in the file ``fodo.lat``.  We first define the relevant constants::
  
   focus    :=   7                   ! [m]     : focal length of equivalent, thin quad
   sepn     :=  10                   ! [m]     : distance between quad centers
   length   :=   2.0                  ! [m]     : quadrupole length
   strength := 1/(focus*length)       ! [m**-2] : quadrupole strength
                                     !         :   = B'/brho

Next, we use these constants to define the elements of the lattice::

   o: drift, l=( sepn - length )
   f: quadrupole, l=length, k1=strength
   d: quadrupole, l=length, k1=(-strength)

Such elements are stitched together to form the lattice,  named appropriately ``fodo``, in one line::

   fodo:  line=( f, o, d, o )

Multiple such latices can be assembled together to make a more complete and complex machine.  We are not quite done yet:
Synergia and CHEF are more than just an optic code: Processes involving real EM fields will be simulated. Non-relativistic
effects could be considered. So a particle and an energy for the reference particle (i.e., the particle staying at the
center of the center of the 6D phase space) needs to be defined::

   beam, particle=proton, energy=1.5
   
This unique line in the lattice file can be inserted anywhere in the file, it makes sense to place at the beginning or at
the end, for sake of clarity. We are now ready to use this lattice. But, prior to run something, details on what we plan to
simulate must be defined.  This is done in two separate user defined Python scripts: the ``fodo_option.py`` and ``fodo.py``.
The first one can be viewed as a simple list of good-old fashion data statements specifying the details about the
simulation.  The second one is a list of action that takes while simulating the problem. To avoid unnecessary verbosity in
this document, let us refrain stating what all these options are, and briefly comment on the way they are defined. Theses
option are defined via a ``synergia_workflow``::

   opts = synergia_workflow.Options("fodo")
   opts.add("radius", 0.1, "aperture radius [m]", float)

A simple Python container class has the name in instantiated.  The option named ``radius`` has the value 0.1 meters, it's
brief definition appears in the thrid argument and it's type is single precision (32 bit, on most machines) floating point. 
It implies that any particles found at a radial distance of more than 10 cm will be lost. 

These options are then used in the main script in file ``fodo.py``.  Our simulation plan starts by setting this aperture
restriction::
    
   for elem in lattice.get_elements():
       elem.set_double_attribute("aperture_radius", opts.radius)


         

A. Via The C++ Interface. 
^^^^^^^^^^^^^^^^^^^^^^^^^


.. _Syphers: http://home.fnal.gov/~syphers/Education/uspas/USPAS08/
