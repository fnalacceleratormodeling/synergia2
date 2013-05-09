Tutorial
========

A first example: simple_fodo_1
------------------------------

In this example, we will simulate a matched bunch passing through a FODO
cell four times. The FODO cell is defined by the file :file:`fodo.lat`, 
which uses Mad8 format.

:file:`fodo.lat`:

.. literalinclude:: ../../examples/fodo_simple1/fodo.lat

The simulation itself is defined by the Python script 
:file:`fodo_simple1.py`.

:file:`fodo_simple1.py`:

.. literalinclude:: ../../examples/fodo_simple1/fodo_simple1.py

The :cpp:class:`Mad8Reader` class reads lattice files in Mad8 format
and produces an object of type :cpp:class:`Lattice`.