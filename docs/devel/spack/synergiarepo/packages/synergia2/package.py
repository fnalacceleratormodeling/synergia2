from spack import *

class Synergia2(CMakePackage):
    """Synergia: an accelerator simulation framework"""
    homepage = "https://cdcvs.fnal.gov/redmine/projects/synergia2"
    url      = "http://compacc.fnal.gov/~amundson/synergia2-2016.03.01.00.tar.bz2"
    
    version('devel', git='http://cdcvs.fnal.gov/projects/synergia2',
            branch='devel')
    
    depends_on("cmake")
    depends_on("mpi", type=('build', 'run'))
    depends_on("chef@devel", type=('build', 'run'))
    depends_on("boost@:1.62+filesystem+iostreams+python+regex+serialization+system+test",
            type=('build', 'run'))
    depends_on("fftw+mpi", type=('build', 'run'))
    depends_on("gsl", type=('build', 'run'))
    depends_on("hdf5+cxx", type=('build', 'run'))
    depends_on("eigen~metis~mpfr~scotch~suitesparse", type=('build', 'run'))
    depends_on("python", type=('build', 'run'))
    depends_on("py-mpi4py", type=('build', 'run'))
    depends_on("py-nose", type=('build', 'run'))
    depends_on("py-pyparsing", type=('build', 'run'))
    depends_on("py-numpy", type=('build', 'run'))

    def setup_environment(self, spack_env, run_env):
#        site_packages_subdir = join_path( ('python' + 
#            str(self.spec['python'].version.up_to(2))),
#            'site-packages')
#        site_packages_dir = join_path(self.spec.prefix.lib,
#                site_packages_subdir)
        site_packages_dir = self.prefix.lib
        spack_env.prepend_path('PYTHONPATH', site_packages_dir)
        run_env.prepend_path('PYTHONPATH', site_packages_dir)

