Using an internal lattice definition: fodo_simple2
==================================================

This example differs from :ref:`section-fodo_simple1` only in the way
the lattice is defined. Instead of reading from a lattice file, we define
the lattice file directly in Synergia.

The simulation itself is defined by the Python script 
:file:`fodo_simple2.py`.

:file:`fodo_simple2.py`:

.. literalinclude:: ../../examples/fodo_simple2/fodo_simple2.py

Lattice element definitions
---------------------------

The Synergia :cpp:class:`Lattice_element` class contains a general set of
information about an element. Each element has a name, a type, a set of
attributes with numerical values and a set of attributes with string values.
Here we define focusing and defocusing quadrupoles (:code:`f` and :code:`d`,
respectively) and a drift (:code:`o`) with length and strength parameters
as used in Mad8.

Lattice definition
------------------

The lattice element parameters are only given meaning when they are converted
into an internal implementation representation through an :cpp:class:`Element_adaptor`.
When we define the :cpp:class:`Lattice`, we give a name for the lattice and
an object of type :cpp:class:`Element_adaptor_map`, which contains a mapping from
element types to objects of type :cpp:class:`Element_adaptor`. Here we tell the lattice
that the elements are to be interpreted as Mad8 elements by passing it an object of type
:cpp:class:`Mad8_adaptor_map`. There are built-in
adaptor maps for Mad8 and MadX elements. The adaptor maps can be extended to accommodate
new element types, or modified to change the implementation behavior of elements.

The :cpp:class:`Lattice` class contains an ordered list of elements. Each element is unique. 
The uniqueness is enforced by storing a *copy* a copy of an element when it is appended 
to the lattice. This means that modifying, say, the :code:`f` object *after* it is appended
to the lattice will not affect the :code:`f` in the lattice itself. It also means that
the lattice contains two different drifts named "o". (In practice, it would have been better
to name the first drift "o1" and the second drift "o2".

Reference particle definition
-----------------------------

In :file:`fodo.lat` the beam energy was defined by the :code:`BEAM` statement.
Within Synergia, the :cpp:class:`Lattice` class contains an object of type 
:cpp:class:`Reference_particle`. Dealing with the various relations between relativistic
velocity, momentum and energy is made simpler by the use of the :cpp:class:`Four_momentum`
class, which does all the relevant conversions.

Remaining definitions and execution
-----------------------------------

Now that the lattice object has been properly defined, the remainder of the example
is exactly the same as :ref:`section-fodo_simple1`.
 
