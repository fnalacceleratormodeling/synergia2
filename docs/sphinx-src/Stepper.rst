Stepper Objects
===============

A step can take place over a fixed distance or a given element of the machine. The
number of steps is always declared at the onset of the simulation.  While the fixed longitudinal distance 
criteria is  more appropriate for simple and highly regular lattices, having the possibility to cover a specific
element in a given step allows for a more precise simulation of the various
collective effects that can take place in a specified environment. 

Default Steppers object are described in the appropriate classes listed in the simulation_ section. 
The user can also define his own stepper, should specialized process be simulated for a given type of
elements. 

.. _simulation: ./simulation.html

