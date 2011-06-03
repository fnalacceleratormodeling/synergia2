Stepper Objects
===============

A step can take place over a fixed distance or a given element of the machine. In the former case, the
number of steps per turn in the machine is decalred at the onset of the simulation. In the latter case,
Synergia will go one element at a time, allowing for a more precise simulation of the various
collective effects that can take place in a specified environment. 

Default Steppers object are described in the appropriate clases listed in the simulation_ section. 
The user can also define his own stepper, should specialized process be simulated for a given type of
elements. 

.. _simulation: ./simulation.html

