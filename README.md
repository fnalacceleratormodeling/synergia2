
# Build Instructions

## 0. Path for non-system packages

    export LOCAL_ROOT=/path/to/local/packages


## 1. General Linux:

    CC=gcc CXX=g++ \
    cmake -DCMAKE_INSTALL_PREFIX=$LOCAL_ROOT \
      -DFFTW3_LIBRARY_DIRS=$LOCAL_ROOT/lib \
      -DEIGEN3_INCLUDE_DIR=$LOCAL_ROOT/include/eigen3 \
      -DHDF5_ROOT=$LOCAL_ROOT \
      -DCMAKE_BUILD_TYPE=Release \
      -DKokkos_ENABLE_OPENMP=on \
      /path/to/synergia/

Kokkos options:

    cmake -DKokkos_ENABLE_OPENMP=on
    cmake -DKokkos_ENABLE_CUDA=off

Enable Python bindings:

    cmake -DBUILD_PYTHON_BINDINGS=on

Enable simple timer profiling:

    cmake -DSIMPLE_TIMER=on


## 2. Cori - KNL:

    module load cmake
    module switch craype-haswell craype-mic-knl
    module load cray-fftw
    module load python
    module load gsl
    module load cray-hdf5
    module unload craype-hugepage2M

    export CRAYPE_LINK_TYPE=dynamic

    CC=cc CXX=CC cmake -DEIGEN3_INCLUDE_DIR:PATH=~/local/include/eigen3 -DFFTW3_LIBRARY_DIRS:PATH=${FFTW_ROOT}/lib -DKokkos_ENABLE_OPENMP=on ../synergia2/


## 3. Power9:

    export PATH=$LOCAL_ROOT/bin:/usr/local/gnu8/gcc-8.3.0/bin:/usr/local/openmpi-4.0.2-gcc-8.3.0/bin:/usr/local/cmake3/3.16.2/bin:/usr/local/cuda/bin:$PATH
    export LD_LIBRARY_PATH=$LOCAL_ROOT/lib:$LOCAL_ROOT/lib64:/usr/local/gnu8/gcc-8.3.0/lib:/usr/local/gnu8/gcc-8.3.0/lib64:/usr/local/cuda/lib64:$LD_LIBRARY_PATH

    CC=/usr/local/gnu8/gcc-8.3.0/bin/gcc CXX=/usr/local/gnu8/gcc-8.3.0/bin/g++ \
    cmake -DEIGEN3_INCLUDE_DIR=/data/qlu/local/include/eigen3 \
      -DFFTW3_LIBRARY_DIRS=/data/qlu/local/lib \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_PYTHON_BINDINGS=off \
      -DKokkos_ENABLE_OPENMP=off \
      -DKokkos_ENABLE_CUDA=on \
      -DCMAKE_CXX_COMPILER=nvcc_wrapper \
      -DCMAKE_CXX_FLAGS="-arch=sm_70" 
      ../synergia2/


Power9 enable Python:

    cmake -DBUILD_PYTHON_BINDINGS=on -DPYTHON_EXECUTABLE=$LOCAL_ROOT/bin/python3


Build python3 on Power9 (libffi needd to enable _ctypes):

    wget python.tar.gz && tar xf
    cd python
    ./configure --prefix=$LOCAL_ROOT
    CFLAGS='-I$LOCAL_ROOT/include' LDFLAGS='-L$LOCAL_ROOT/lib64 -lffi' make -j 32
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



## 4. MacOS:

    brew install gcc hdf5 eigen fftw3

MacOS with gcc:

    CC=gcc-9 CXX=g++-9 cmake -DKokkos_ENABLE_OPENMP=on -DCMAKE_OSX_SYSROOT="/" -DCMAKE_OSX_DEPLOYMENT_TARGET="" /path/to/synergia/

MacOS with apple clang (openmp missing, possible with '-Xpreprocessor -fopenmp -lomp'):

    cmake /path/to/synergia

4.3 MacOS with Python 3 (homebrew):

    cmake -DPYTHON_EXECUTABLE=/usr/local/bin/python3
