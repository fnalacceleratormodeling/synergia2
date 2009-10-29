project(SYNERGIA2)
cmake_minimum_required(VERSION 2.6)
enable_testing()

# python
INCLUDE(FindPythonInterp)
INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})

# boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD ON)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework)

#~ message("Boost_INCLUDE_DIRS = '${Boost_INCLUDE_DIRS}'")
#~ message("Boost_LIBRARIES = '${Boost_LIBRARIES}'")
#~ message("Boost_LIBRARY_DIRS = '${Boost_LIBRARY_DIRS}'")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# chef
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindCHEF.cmake)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

# numpy
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindNumpy.cmake)
include_directories(${NUMPY_INCLUDE_DIR})

# mpi4py
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindMpi4py.cmake)
include_directories(${MPI4PY_INCLUDE_DIR})

# misc
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindFFTW2.cmake)

# for eclipse only
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia-local.in" 
    "${SYNERGIA2_BINARY_DIR}/synergia-local" IMMEDIATE)

include_directories(${SYNERGIA2_SOURCE_DIR})

add_custom_target(doxy)
add_custom_target(sphinxdoc)
add_dependencies(sphinxdoc doxy)

add_subdirectory(array_nd)
add_subdirectory(utils)
link_directories(utils)
add_subdirectory(components)
add_subdirectory(synergia)
add_subdirectory(s2_fish)
add_subdirectory(impact)

add_subdirectory(octapy)

add_subdirectory(doc)