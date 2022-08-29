import sys
sys.path.append('..')
sys.path.append('../../foundation')
sys.path.append('../../convertors')

from mpi4py import MPI
from foundation import Reference_particle, Four_momentum
from bunch import Bunch
from bunch import Diagnostics_writer, Multi_diagnostics_writer

def test_construct():
    multi_diagnostics_writer = Multi_diagnostics_writer()

def test_append():
    multi_diagnostics_writer = Multi_diagnostics_writer()
    diagnostics_writer1 = Diagnostics_writer()
    multi_diagnostics_writer.append(diagnostics_writer1)
    diagnostics_writer2 = Diagnostics_writer()
    multi_diagnostics_writer.append(diagnostics_writer2)
