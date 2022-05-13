# Build Instructions

## 0a. Path for non-system packages

    export LOCAL_ROOT=/path/to/local/packages

## 0b. Clone repository

Git 2.13 and later,

    git clone -b devel3 --recurse-submodules https://github.com/fnalacceleratormodeling/synergia2.git

Git 1.6 and later,

    git clone -b devel3 --recursive https://github.com/fnalacceleratormodeling/synergia2.git

Or manually init and update submodules,

    git clone -b devel3 https://github.com/fnalacceleratormodeling/synergia2.git
    cd synergia2
    git submodule update --init --recursive

## 0c. Update local repo

In one command

    git pull --recurse-submodules

or,

    git pull
    git submodule update

## 1. General Linux:

    CC=gcc CXX=g++ \
    cmake -DCMAKE_INSTALL_PREFIX=$LOCAL_ROOT \
      -DCMAKE_BUILD_TYPE=Release \
      -DFFTW3_LIBRARY_DIRS=$LOCAL_ROOT/lib \
      -DHDF5_ROOT=$LOCAL_ROOT \
      -DKokkos_ENABLE_OPENMP=on \
      /path/to/synergia/

### Available Build Options

Kokkos options:

    cmake -DKokkos_ENABLE_OPENMP=on/off
    cmake -DKokkos_ENABLE_CUDA=on/off
    ...

Synergia2 requires an OpenMP-capable C++ compiler (when compiling for host-only platforms). On MacOS, AppleClang is not OpenMP-enabled by default and one has to install libomp via brew. Note that we carry [a bug fix](https://github.com/fnalacceleratormodeling/kokkos/commit/b2057c3bd455d3d9af8993187f5add348fff1549) for this case in the vendored copy of Kokkos. You can choose not to use it by passing `-DUSE_EXTERNAL_KOKKOS=on` and provide an external installation of Kokkos.

Enable/disable Python bindings:

    cmake -DBUILD_PYTHON_BINDINGS=on|off  # default is on

Vectorization flags (on M1 Mac must set `-DGSV=DOUBLE`):

    cmake -DGSV=DOUBLE|SSE|AVX|AVX512

Enable/disable simple timer profiling:

    cmake -DSIMPLE_TIMER=on|off # default is off

### Option to enable OpenMP backend for kokkos (host-only build)

    cmake -DENABLE_OPENMP=on|off # default is off

### Options to enable CUDA backend for kokkos

    cmake -DENABLE_CUDA=on|off

Please ensure that both OpenMP and CUDA are not enabled together.

Paddings need to be turned off in the CUDA build due to a Kokkos bug https://github.com/kokkos/kokkos/issues/2995

    cmake -DALLOW_PADDING=on|off

### Other dependencies

We allow external installations of cereal and pybind11 (for building the python bindings) if preferred by the user. These can be passed by the flags `-DUSE_EXTERNAL_CEREAL/PYBIND11=on` and ensuring that CMake can find them. The defaul behavior is to fetch the latest versions of these packages during configuration.

When running on a cluster that uses slurm as the resource allocation manger/job launch tool, we assume that the default launch flags for MPI programs is `srun --mpi=pmix_v3`. If you are using a different version of pmi2/pmix, please pass `-DSRUN_MPI_PMIX=...`.

## 2. Ubuntu 20.04 LTS

Out general philosophy is to use the package manager to install as many
requirements as possible.

    sudo apt install \
        cmake \
        cython3 \
        g++ \
        hdf5-tools \
        libfftw3-dev \
        libfftw3-mpi-dev \
        libgsl-dev \
        libhdf5-dev \
        libhdf5-openmpi-dev \
        libopenmpi-dev \
        libpython3-dev \
        pkg-config \
        python3-dev \
        python3-matplotlib \
        python3-mpi4py \
        python3-pkgconfig \
        python3-pyparsing \
        python3-pytest \
        python3-tables \
        python3.8-venv


    # We recommend a python virtual environment for managing module versions.
    python3 -m venv --system-site-packages synergia-env
    source synergia-env/bin/activate

    # We do not install h5py using pip because the pip-installed version may not
    # not use a version of HDF5 that matches what we have from the package manager.
    i# Instead, we build our own from source.
    mkdir tmp
    cd tmp
    wget https://github.com/h5py/h5py/archive/refs/tags/3.6.0.tar.gz # check for the most recent version
    tar xf 3.6.0.tar.gz
    cd h5py-3.6.0
    HDF5_DIR=/usr/local python setup.py install
    cd ../..
    rm -r tmp/

    # We build using the system compilers which are recent enough for our needs
    cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on /path/to/synergia/

## 3. macOS (Intel and M1 Mac)

Our general philosophy is to use Homebrew to install as many requirements as possible.
See [https://brew.sh](https://brew.sh) for instructions on the installation and use of Homebrew.

    # Note that the Homebrew installation of hdf5 does not (at the time of this writing)
    # support MPI parallelism.
    brew install gsl hdf5 fftw libomp ninja

    # We recommend a python virtual environment for managing module versions.
    python3 -m venv synergia-env
    source synergia-env/bin/activate
    python3 -m pip install --upgrade pip
    python3 -m pip install cython matplotlib mpi4py numpy pyparsing pytest

    # We do not install h5py using pip because the pip-installed version may not
    # not use a version of HDF5 that matches what we have from the package manager.
    # Instead, we build our own from source.
    mkdir tmp
    cd tmp
    wget https://github.com/h5py/h5py/archive/refs/tags/3.6.0.tar.gz # check for the most recent version
    tar xf 3.6.0.tar.gz
    cd h5py-3.6.0
    HDF5_DIR=/usr/local python setup.py install
    cd ../..
    rm -r tmp/

### macOS with apple clang

    # We do not recommend using /usr/local/ as your installation target. While this is the default, this will mix your Synergia installation
    # with the tools installed using Homebrew -- but Homebrew will not know how to update Synergia.
    cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on -DCMAKE_INSTALL_PREFIX=/path/to/install/target /path/to/synergia/
    

Note that updates of Homebrew-installed packages can invalidate an existing build.
Often rerunning the build (with `make` or `ninja`) will be sufficient.
Sometimes, however, this results in re-building or even re-configuration failures.
The solution in this case is a complete clean re-installation.
Delete all files (including hidden files) in the build directory, and re-run `cmake`.

### macOS with gcc

We do not recommend the use of the the GNU compiler suite to build Synergia on macOS.
This can lead to incompatibilities between C++ libraries that are part of the OS or part of Homebrew on one hand, and C++ libraries built with `g++` on the other.

    # The current version of the Homebrew GCC formula at the time of this writing installs
    # g++-11.
    CC=gcc-11 CXX=g++-11 cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on -DBUILD_PYTHON_BINDINGS=on /path/to/synergia/

# Obsolete instructions

The following platforms have been supported in the past, but their support has lapsed.
These instructions may no longer work.
If you are trying to install on one of these platforms and encounter a problem, please contact us by creating an issue on GitHub, to start a discussion on whether reviving support is feasible.

## 4. Cori - KNL (obsolete, needs update):

    module load cmake
    module switch craype-haswell craype-mic-knl
    module load cray-fftw
    module load python
    module load gsl
    module load cray-hdf5
    module unload craype-hugepage2M

    export CRAYPE_LINK_TYPE=dynamic

    CC=cc CXX=CC cmake -DEIGEN3_INCLUDE_DIR:PATH=~/local/include/eigen3 -DFFTW3_LIBRARY_DIRS:PATH=${FFTW_ROOT}/lib -DKokkos_ENABLE_OPENMP=on ../synergia2/


## 5. Power9 (obsolete, needs update):

    export SYN_SRC=/path/to/synergia
    export LOCAL_ROOT=/data/qlu/local

    export PATH=$LOCAL_ROOT/bin:/usr/local/gnu8/gcc-8.3.0/bin:/usr/local/openmpi-4.0.2-gcc-8.3.0/bin:/usr/local/cmake3/3.16.2/bin:/usr/local/cuda-11.1/bin:$PATH
    export LD_LIBRARY_PATH=$LOCAL_ROOT/lib:$LOCAL_ROOT/lib64:/usr/local/gnu8/gcc-8.3.0/lib:/usr/local/gnu8/gcc-8.3.0/lib64:/usr/local/cuda-11.1/compat:$LD_LIBRARY_PATH

    CXX=/usr/local/gnu8/gcc-8.3.0/bin/g++ \
    cmake -DEIGEN3_INCLUDE_DIR=$LOCAL_ROOT/include/eigen3 \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_PYTHON_BINDINGS=off \
      -DKokkos_ENABLE_OPENMP=off \
      -DKokkos_ENABLE_CUDA=on \
      -DCMAKE_CXX_COMPILER=$SYN_SRC/src/synergia/utils/kokkos/bin/nvcc_wrapper \
      -DCMAKE_CXX_FLAGS="-arch=sm_70"
      ../synergia2/


Power9 enable Python:

    cmake -DBUILD_PYTHON_BINDINGS=on -DPYTHON_EXECUTABLE=$LOCAL_ROOT/bin/python3


Build python3 on Power9 (libffi needd to enable _ctypes):

    wget python.tar.gz && tar xf
    cd python
    LDFLAGS='-L${LOCAL_ROOT}/lib64' CFLAGS='-I${LOCAL_ROOT}/include' ./configure --prefix=$LOCAL_ROOT  --enable-optimizations --enable-shared
    make -j 32
    make install

    wget mpi4py.tar.gz && tar xf
    cd mpi4py
    python3 setup.py build
    python3 setup.py install

    wget Cython.tar.gz && tar xf
    cd Cython
    python3 setup.py install

    git clone https://github.com/numpy/numpy.git numpy
    cd numpy
    python3 setup.py install

