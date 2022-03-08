import pytest


def test_synergia_version():
    """Make sure that the version number generated at configuration time
    is the expected one.
    """
    from synergia import version

    assert version.__version__ == "2018.02.20.00"


def test_python_version():
    """Make sure that the Python interpreter version as found by CMake at
    configuration time matches the version found by the runtime, when this
    test is executed.
    """
    from synergia import version
    from sys import version_info as vi

    assert version.python_interp_version == f"{vi.major}.{vi.minor}.{vi.micro}"


def test_mpi_version():
    """Make sure that the MPI library found at runtime is the same as that
    which was found by CMake at configuration time.
    """
    from synergia import version
    from mpi4py import MPI

    # For some reason, the module function MPI.Get_library_version returns
    # a string will a trailing null character, which we remove.
    assert version.mpi_library_version == MPI.Get_library_version().rstrip("\0")

def test_hdf5_version():
    from synergia import version
    from h5py import version as v
    assert version.hdf5_library_version == v.hdf5_version
    