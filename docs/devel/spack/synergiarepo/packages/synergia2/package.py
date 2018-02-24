from spack import *

class Synergia2(CMakePackage):
    """Synergia: an accelerator simulation framework"""
    homepage = "https://cdcvs.fnal.gov/redmine/projects/synergia2"
    url      = "http://compacc.fnal.gov/~amundson/synergia2-2016.03.01.00.tar.bz2"
    
    version('devel', git='http://cdcvs.fnal.gov/projects/synergia2',
            branch='devel')
    
    depends_on("cmake")
    depends_on("mpi")
    depends_on("chef@devel")
    depends_on("boost@:1.62+filesystem+iostreams+python+regex+serialization+system+test")
    depends_on("fftw+mpi")
    depends_on("gsl")
    depends_on("hdf5+cxx")
    depends_on("eigen~metis~mpfr~scotch~suitesparse")
    depends_on("python")
    depends_on("py-mpi4py")
    depends_on("py-nose")
    depends_on("py-pyparsing")
    depends_on("py-numpy")

