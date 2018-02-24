from spack import *

class Chef(CMakePackage):
    """CHEF"""
    homepage = "https://cdcvs.fnal.gov/redmine/projects/chef"

    version('devel', 
            git='http://cdcvs.fnal.gov/projects/accelerator-modeling-chef',
            branch='devel-synergia')

    depends_on("cmake")
    depends_on("boost@:1.62+filesystem+iostreams+python+regex+serialization+system+test")
    depends_on("fftw+mpi")
    depends_on("python")
    depends_on("py-numpy")

    def cmake_args(self):
        spec = self.spec
        args = ['-DBUILD_PARSER_MODULES=OFF', 
                '-DBUILD_PYTHON_BINDINGS=ON',
                '-DBUILD_SHARED_LIBS=ON']
        return args

    def setup_environment(self, spack_env, run_env):
        site_packages_subdir = join_path( ('python' + 
            str(self.spec['python'].version.up_to(2))),
            'site-packages')
        site_packages_dir = join_path(self.spec.prefix.lib,
                site_packages_subdir)
        spack_env.prepend_path('PYTHONPATH', site_packages_dir)
        run_env.prepend_path('PYTHONPATH', site_packages_dir)
