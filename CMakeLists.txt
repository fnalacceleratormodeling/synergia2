if("${CMAKE_CUDA_COMPILER}" MATCHES "clang")
  cmake_minimum_required(VERSION = 3.18)
else()
  cmake_minimum_required(VERSION 3.16)
endif()

project(SYNERGIA3 LANGUAGES CXX C)

if(${CMAKE_VERSION} VERSION_GREATER "3.22.0")
  cmake_policy(SET CMP0127 OLD) # remove this when we move to CMake 3.22
endif()

# Syenrgia3 Version
set(SYNERGIA3_VERSION_YEAR 2024)
set(SYNERGIA3_VERSION_MONTH 04)
set(SYNERGIA3_VERSION_DAY 17)
if(DEFINED ENV{GITHUB_SHA})
  set(SYNERGIA3_GIT_HASH_FULL $ENV{GITHUB_SHA})
  string(SUBSTRING ${SYNERGIA3_GIT_HASH_FULL} 0 7 SYNERGIA3_GIT_HASH)
else()
  # Get the latest abbreviated commit hash of the working branch
  # https://jonathanhamberg.com/post/cmake-embedding-git-hash/
  execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
    OUTPUT_VARIABLE SYNERGIA3_GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(SYNERGIA3_VERSION
    "${SYNERGIA3_VERSION_YEAR}.${SYNERGIA3_VERSION_MONTH}.${SYNERGIA3_VERSION_DAY}-${SYNERGIA3_GIT_HASH}"
)
message(STATUS "Preparing to install Synergia3 version: ${SYNERGIA3_VERSION}")

# generate compile_commands.json
if(${CMAKE_GENERATOR} STREQUAL "Unix Makefiles" OR ${CMAKE_GENERATOR} STREQUAL
                                                   "Ninja")
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# project specific cmake modules
set(CMAKE_MODULE_PATH "${SYNERGIA3_SOURCE_DIR}/CMake")
include(CopyFile)
include(MPITest)
include(AddPythonTest)

# Make sure we know which Kokkos back end we're using. If none is specified on
# the command line, the default is OpenMP.
if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "")
  set(ENABLE_KOKKOS_BACKEND OpenMP)
  message(STATUS "OpenMP backend for Kokkos is enabled.")
endif()

# Make sure the ENABLE_KOKKOS_BACKEND is set to one of the allowed values, and
# nothing else.
if((NOT "${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
   AND (NOT "${ENABLE_KOKKOS_BACKEND}" STREQUAL "OpenMP"))
  message(FATAL_ERROR "ENABLE_KOKKOS_BACKEND must be one of CUDA or OpenMP")
else()
  message(STATUS "ENABLE_KOKKOS_BACKEND is ${ENABLE_KOKKOS_BACKEND}")
endif()

# Set preprocessor flags we require
if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "OpenMP")
  message(STATUS "Adding SYNERGIA_ENABLE_OPENMP to compiler definiions")
  add_compile_definitions(SYNERGIA_ENABLE_OPENMP=1)
  set(Kokkos_ENABLE_CUDA
      Off
      CACHE BOOL "If true-ish use CUDA backend for Kokkos" FORCE)
  set(Kokkos_ENABLE_OPENMP
      On
      CACHE BOOL "If true-ish use OpenMP backend for Kokkos" FORCE)
endif()

if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
  message(STATUS "Adding SYNERGIA_ENABLE_CUDA to compiler definiions")
  add_compile_definitions(SYNERGIA_ENABLE_CUDA=1)
  set(Kokkos_ENABLE_CUDA
      On
      CACHE BOOL "If true-ish use CUDA backend for Kokkos" FORCE)
  set(Kokkos_ENABLE_OPENMP
      Off
      CACHE BOOL "If true-ish use OpenMP backend for Kokkos" FORCE)
endif()

# Set up to assure fully-linked libraries.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(COMPILER_SPECIFIC_CXX_FLAGS
      "-Wno-#pragma-messages -Wno-potentially-evaluated-expression -Wno-register"
  )
  message(STATUS "Using clang++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Apples's clang
  set(COMPILER_SPECIFIC_CXX_FLAGS
      "-Wno-#pragma-messages -Wno-potentially-evaluated-expression -Wno-register"
  )
  message(STATUS "Using Apple clang++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-error=unused-result -Wno-register")
  if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(LINKER_OPTIONS "LINKER:-z,defs")
  endif()
  message(STATUS "Using g++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-error=unused-result")
  set(LINKER_OPTIONS "LINKER:-z,defs")
  message(STATUS "Using classic Intel C++")
else()
  message(
    WARNING
      "I can not tell what compiler is being used. Will proceed with generic settings."
  )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMPILER_SPECIFIC_CXX_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS
    "${CMAKE_SHARED_LINKER_FLAGS} ${COMPILER_SPECIFIC_SHARED_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS
    "${CMAKE_MODULE_LINKER_FLAGS} ${COMPILER_SPECIFIC_SHARED_LINKER_FLAGS}")

# options
option(SIMPLE_TIMER "use the simple timer 2 profiling" OFF)
option(BUILD_PYTHON_BINDINGS "build the Python bindings of Synergia" ON)
option(ALLOW_PADDING "enable particle array padding for memory alignment" OFF)
option(USE_OPENPMD_IO
       "enable the use of openPMD for parallel I/O, requires parallel HDF5" OFF)
option(USE_EXTERNAL_KOKKOS "Use external build of kokkos" OFF)
option(USE_EXTERNAL_CEREAL "Use external build of cereal" OFF)
option(
  USE_EXTERNAL_PYBIND11
  "Use external build of pybind11, only applicable when building python bindings"
  OFF)
option(USE_EXTERNAL_OPENPMD "Use external build of openPMD" OFF)
option(BUILD_EXAMPLES "build examples" ON)
option(
  BUILD_FD_SPACE_CHARGE_SOLVER
  "build the finite difference based space charge solver, requires PETSc & PkgConfig"
  OFF)

# Convert the variables to uppercase to ensure generator expressions work as
# expected
string(TOUPPER ${USE_OPENPMD_IO} USE_OPENPMD_IO)
string(TOUPPER ${BUILD_FD_SPACE_CHARGE_SOLVER} BUILD_FD_SPACE_CHARGE_SOLVER)

# PETSc
if(${BUILD_FD_SPACE_CHARGE_SOLVER})
  add_definitions(-DBUILD_FD_SPACE_CHARGE_SOLVER)
  # set root of location to find PETSc's pkg-config
  set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
  set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH})
  find_package(PkgConfig REQUIRED)
  pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc>=3.19.0)
  if(${PETSC_VERSION} STRGREATER_EQUAL 3.20.0)
    message(STATUS "Found PETSc version ${PETSC_VERSION}")
  else()
    message(
      STATUS
        "Found PETSc version ${PETSC_VERSION}, using default CPU based LU facorization, upgrade to version 3.20.0 and install STRUMPACK for enhanced performance on GPUs when using fixed domains"
    )
  endif()
endif()

# srun --mpi=SRUN_MPI_PMIX
set(SRUN_MPI_PMIX
    "pmix_v3"
    CACHE STRING "PMI2/PMIx version to use with srun")

# only build the shared lib
option(BUILD_SHARED_LIBS "Build shared libraries" TRUE)

# option for gsvector setting
set(GSV
    "AVX"
    CACHE STRING "set gsvector vectorization type (DOUBLE|SSE|AVX|AVX2|AVX512)")

# find MPI -- do not build the C++ bindings
set(MPI_CXX_SKIP_MPICXX
    TRUE
    CACHE BOOL "disable the MPI2-C++ bindings")
set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI REQUIRED COMPONENTS C CXX) # CXX means we need to call C-APIs
                                            # in C++ code
add_compile_definitions(${MPI_CXX_COMPILE_DEFINITIONS})
message("-- MPI_VERSION is: ${MPI_VERSION}")

# launch tests with additional --mpi=${SRUN_MPI_PMIX} flag!
if(${MPIEXEC_EXECUTABLE} MATCHES "srun")
  set(MPIEXEC_EXECUTABLE ${MPIEXEC_EXECUTABLE} "--mpi=${SRUN_MPI_PMIX}")
endif()

# We require OpenMP support
find_package(OpenMP REQUIRED)

# Follow official kokkos guidelines
# https://github.com/kokkos/kokkos/wiki/Compiling#42-using-general-cmake-build-system
if(USE_EXTERNAL_KOKKOS)
  find_package(Kokkos REQUIRED)
  message(STATUS "Using external installation of kokkos")
  set(kokkos_libs "Kokkos::kokkos")
  if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
    project(SYNERGIA3 LANGUAGES CXX C CUDA)
  endif()
else()
  message(STATUS "Using vendored copy of kokkos")
  if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
    set(Kokkos_ENABLE_CUDA_LAMBDA
        ON
        CACHE BOOL "enable CUDA LAMBDA feature" FORCE)
    message(STATUS "Using vendored copy of kokkos with CUDA")
  endif()
  add_subdirectory(src/synergia/utils/kokkos)
  include_directories(${Kokkos_INCLUDE_DIRS_RET})
  set(kokkos_libs kokkos)
endif()

# Check version of kokkos!
message("-- Using Kokkos version: ${Kokkos_VERSION}")
if(Kokkos_VERSION LESS 4.0.00)
  message(FATAL_ERROR "Need at least Kokkos 4.0.00!")
endif()

# Catch2-v3 for testing! This has to happen after $kokkos_libs is defined!
if(USE_EXTERNAL_CATCH2)
  find_package(Catch2 3 REQUIRED CONFIG)
  message(STATUS "Using external catch2-v3")
else()
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4)
  FetchContent_MakeAvailable(Catch2)
endif()
set(testing_libs
    Catch2::Catch2 ${kokkos_libs} MPI::MPI_C
    $<$<STREQUAL:${BUILD_FD_SPACE_CHARGE_SOLVER},ON>:PkgConfig::PETSC>)
set(test_main "${CMAKE_SOURCE_DIR}/src/synergia/utils/catch_test_main.cc")
enable_testing()

# HDF5
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5 REQUIRED)

if(HDF5_IS_PARALLEL)
  message(STATUS "Using parallel hdf5")
  add_definitions(-DUSE_PARALLEL_HDF5)
  if(USE_OPENPMD_IO)
    message(STATUS "Building openPMD I/O capabilities")
    # OpenPMD
    if(USE_EXTERNAL_OPENPMD)
      find_package(openPMD CONFIG REQUIRED)
      message(STATUS "using external installation of openPMD")
    else()
      message(STATUS "fetching OpenPMD")
      include(FetchContent)
      set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
      set(openPMD_BUILD_CLI_TOOLS OFF)
      set(openPMD_BUILD_EXAMPLES OFF)
      set(openPMD_BUILD_SHARED_LIBS OFF)
      set(openPMD_BUILD_TESTING OFF)
      set(openPMD_INSTALL OFF)
      set(openPMD_USE_PYTHON ON)
      FetchContent_Declare(
        openPMD
        GIT_REPOSITORY "https://github.com/openPMD/openPMD-api.git"
        GIT_TAG 0.15.2)
      FetchContent_MakeAvailable(openPMD)
    endif()
    add_definitions(-DSYNERGIA_HAVE_OPENPMD)
  else()
    message(STATUS "openPMD I/O capabilties disabled")
  endif()
else(HDF5_IS_PARALLEL)
  message(STATUS "Using serial hdf5")
  message(STATUS "openPMD I/O capabilties disabled")
endif()

# GSL
find_package(GSL REQUIRED)

# FFTW3
if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "OpenMP")
  find_package(FFTW3 REQUIRED)

  if(NOT FFTW3_MPI_FOUND)
    message(FATAL_ERROR "Missing FFTW3-mpi")
  endif()

  if(NOT FFTW3_OMP_FOUND)
    message(FATAL_ERROR "FFTW3 openmp not found")
  else()
    set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_OMP_LIBRARIES}
                                ${FFTW3_LIBRARIES})
  endif()
endif()

# cereal
if(USE_EXTERNAL_CEREAL)
  find_package(cereal REQUIRED)
  message(
    STATUS "Using external cereal with cmake configuration from ${cereal_DIR}")
else()
  message(STATUS "Fetching cereal")
  include(FetchContent)
  set(BUILD_SANDBOX OFF)
  set(SKIP_PERFORMANCE_COMPARISON ON)
  set(SKIP_PORTABILITY_TEST ON)
  set(JUST_INSTALL_CEREAL ON)
  set(CEREAL_INSTALL ON)
  FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    GIT_TAG v1.3.2)
  FetchContent_MakeAvailable(cereal)
endif()

# warnings
string(APPEND ${CMAKE_CXX_FLAGS} "-Wall -Werror")

# enable fast-math, but disable the finite-math-only to avoid the
# std::isfinite() woes
string(APPEND ${CMAKE_CXX_FLAGS} "-ffast-math -fno-finite-math-only")

# for eclipse only (no line wrapping for error messages)
if(CMAKE_COMPILER_IS_GNUCXX)
  string(APPEND ${CMAKE_CXX_FLAGS} "-fmessage-length=0")
endif(CMAKE_COMPILER_IS_GNUCXX)

# python
if(BUILD_PYTHON_BINDINGS)
  message(STATUS "Building python bindings")
  find_package(Python REQUIRED COMPONENTS Interpreter Development)

  # Make sure we have mpi4py available.
  execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import mpi4py"
    RESULT_VARIABLE EXIT_CODE
    OUTPUT_QUIET)

  if(NOT ${EXIT_CODE} EQUAL 0)
    message(
      FATAL_ERROR
        "Python_EXECUTABLE is ${Python_EXECUTABLE}\nThe \"mpi4py\" module is not availale."
    )
  endif()

  # pybind11
  if(USE_EXTERNAL_PYBIND11)
    find_package(pybind11 REQUIRED CONFIG)
    message(STATUS "Using external pybind11")
  else()
    message(STATUS "Fetching pybind11")
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11
      GIT_TAG v2.12.0)
    FetchContent_MakeAvailable(pybind11)
  endif()
else(BUILD_PYTHON_BINDINGS)
  message(STATUS "Do not build python bindings")
endif()

# simple timer
if(SIMPLE_TIMER)
  message(STATUS "Simple timer 2 profiling enabled")
  add_compile_definitions(SIMPLE_TIMER=1)
endif()

# allow padding
if(ALLOW_PADDING)
  message(STATUS "Particle array padding enabled")
else()
  message(STATUS "Particle array padding disabled")
  add_definitions(-DNO_PADDING)
endif()

# GSVector settings
if(GSV STREQUAL "SSE")
  message(STATUS "using SSE for gsvector")
  add_definitions(-DGSV_SSE)
elseif(GSV STREQUAL "AVX")
  message(STATUS "using AVX for gsvector")
  add_definitions(-DGSV_AVX)
elseif(GSV STREQUAL "AVX2")
  message(STATUS "using AVX2 for gsvector")
  add_definitions(-DGSV_AVX2)
elseif(GSV STREQUAL "AVX512")
  message(STATUS "using AVX512 for gsvector")
  add_definitions(-DGSV_AVX512)
else()
  message(STATUS "using double for gsvector")
endif()

# Additional CXXFLAGS
if(EXTRA_CXX_FLAGS)
  string(APPEND CMAKE_CXX_FLAGS ${EXTRA_CXX_FLAGS})
endif()

# External libraries generate many warnings add_definitions(-Wall -pedantic
# -Wno-long-long)

configure_file("${SYNERGIA3_SOURCE_DIR}/synergia-local.in"
               "${SYNERGIA3_BINARY_DIR}/synergia-local" IMMEDIATE)

configure_file("${SYNERGIA3_SOURCE_DIR}/src/local_paths.py.in"
               "${SYNERGIA3_BINARY_DIR}/src/local_paths.py" IMMEDIATE)

include_directories(BEFORE ${SYNERGIA3_SOURCE_DIR}/src)

if(NOT TARGET test-exes)
  add_custom_target(test-exes)
endif()

function(add_test_executable name)
  add_executable(${name} ${ARGN})
  add_dependencies(test-exes ${name})
endfunction()

add_custom_target(
  check
  COMMAND ${CMAKE_CTEST_COMMAND}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Running ctest")
add_dependencies(check test-exes)

set(INCLUDE_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/include/"
    CACHE PATH "include install directory")
set(LIB_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/"
    CACHE PATH "library install directory")
set(BIN_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/bin/"
    CACHE PATH "executable install directory")
set(PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/"
    CACHE PATH "python install directory")

set(SYNERGIA_EXTRA_LIBRARY_PATHS
    ""
    CACHE STRING "additional library paths for synergia executable")
set(SYNERGIA_EXTRA_PYTHON_PATHS
    ""
    CACHE STRING "additional python module paths for synergia executable")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${LIB_INSTALL_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# possibly move the following to using block/endblock when cmake min version
# becomes 3.25 Find Eigen, a dependency we vendor
message(STATUS "Configuring Eigen3 for use in host routines.")
# Eigen tries to determine fortran compiler ID via FC, which can cause issues
# with flang-new, so unset that temporarily.
set(FC_BACKUP $ENV{FC})
unset(ENV{FC})
# https://gitlab.com/libeigen/eigen/-/issues/848
set(eigen_defs
    EIGEN_NO_CUDA
    $<$<STREQUAL:${CMAKE_HOST_SYSTEM_PROCESSOR},ppc64le>:EIGEN_DONT_VECTORIZE>)
add_subdirectory(${SYNERGIA3_SOURCE_DIR}/src/synergia/utils/eigen)
set(ENV{FC} ${FC_BACKUP})
add_subdirectory(src)

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
add_subdirectory(synergia-script-templates)
