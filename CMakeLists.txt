project(SYNERGIA2)

INCLUDE(FindPythonInterp)

INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})

FIND_PACKAGE(Boost)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
set(BOOST_LIBS -lboost_python-gcc)

INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindCHEF.cmake)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindFFTW2.cmake)

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia_local.in" 
    "${SYNERGIA2_BINARY_DIR}/synergia_local" IMMEDIATE)

include_directories(${SYNERGIA2_SOURCE_DIR})

add_subdirectory(array_nd)
add_subdirectory(synergia)
add_subdirectory(s2_fish)
add_subdirectory(impact)

add_subdirectory(octapy)
