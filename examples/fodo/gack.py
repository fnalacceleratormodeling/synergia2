#!/usr/bin/env python3

import sys, os
from mpi4py import MPI
import numpy as np
import synergia


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


def run():

    lattice = get_lattice()
    print('Read lattice, length: ', lattice.get_length(), ', ', len(lattice.get_elements()), ' elements')
    if len(lattice.get_elements()) != 4:
        print('ERROR!!! Wrong number of lattice elements!!!!!')
        for elem in lattice.get_elements():
            elem.print_()

def main():

    print("running fodo.py: my rank =", MPI.COMM_WORLD.Get_rank())
    try:
        run()
    except:
        raise RuntimeError("Failure to launch fodo.run")

main()
