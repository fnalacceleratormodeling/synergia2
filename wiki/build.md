
# Build Instructions

## 0a. Path for non-system packages

    export LOCAL_ROOT=/path/to/local/packages

## 0b. Clone repository

Git 2.13 and later,

    git clone -b devel3 --recurse-submodules https://bitbucket.org/fnalacceleratormodeling/synergia2.git

Git 1.6 and later,

    git clone -b devel3 --recursive https://bitbucket.org/fnalacceleratormodeling/synergia2.git

Or manually init and update submodules,

    git clone -b devel3 https://bitbucket.org/fnalacceleratormodeling/synergia2.git
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

### Avaiable Build Options

Kokkos options:

    cmake -DKokkos_ENABLE_OPENMP=on/off
    cmake -DKokkos_ENABLE_CUDA=on/off
    ...

Enable Python bindings:

    cmake -DBUILD_PYTHON_BINDINGS=on
    
Vectorization flags (on M1 Mac must set `-DGSV=DOUBLE`):

    cmake -DGSV=DOUBLE|SSE|AVX|AVX512

Enable simple timer profiling:

    cmake -DSIMPLE_TIMER=on

### Options for OpenMP only Build

    cmake -DKokkos_ENABLE_OPENMP=on

### Options for GPU/CUDA build

`nvcc` needs to be in path

    export PATH=/usr/local/cuda/bin:$PATH

Kokkos options (it is possible to have both openmp and cuda enabled)

    cmake -DKokkos_ENABLE_OPENMP=on
    cmake -DKokkos_ENABLE_CUDA=on

Use `nvcc_wrapper` as the default cxx compiler, and set the GPU architecture in the `CXX_FLAGS`. `nvcc_wrapper` can be found in Synergia source tree under `src/synergia/utils/kokkos/bin/nvcc_wrapper`

    cmake -DCMAKE_CXX_COMPILER=/path/to/nvcc_wrapper

Paddings need to be turned off in the CUDA build due to a Kokkos bug https://github.com/kokkos/kokkos/issues/2995

    cmake -DALLOW_PADDING=off


## 2. Ubuntu20.04

    sudo apt install libopenmpi-dev libfftw3-mpi-dev libgsl-dev libhdf5-openmpi-dev libpython3-dev
    cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on -DBUILD_PYTHON_BINDINGS=on /path/to/synergia/

## 3. Cori - KNL (obsolete, needs update): 

    module load cmake
    module switch craype-haswell craype-mic-knl
    module load cray-fftw
    module load python
    module load gsl
    module load cray-hdf5
    module unload craype-hugepage2M

    export CRAYPE_LINK_TYPE=dynamic

    CC=cc CXX=CC cmake -DEIGEN3_INCLUDE_DIR:PATH=~/local/include/eigen3 -DFFTW3_LIBRARY_DIRS:PATH=${FFTW_ROOT}/lib -DKokkos_ENABLE_OPENMP=on ../synergia2/


## 4. Power9 (obsolete, needs update):

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



## 5. MacOS (Intel/M1 Mac):

    brew install gcc hdf5 fftw3 libomp
    pip3 install numpy mpi4py pytest pyparsing

MacOS with gcc:

    CC=gcc-9 CXX=g++-9 cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on -DBUILD_PYTHON_BINDINGS=on /path/to/synergia/

MacOS with apple clang:

    cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=on -DBUILD_PYTHON_BINDINGS=on /path/to/synergia/

