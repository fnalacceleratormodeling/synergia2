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

add_subdirectory(s2_fish)

add_python_extension(apply_map apply_map.cc)
target_link_libraries(apply_map  ${CHEF_LIBS} ${BOOST_LIBS})

add_python_extension(chef_propagate chef_propagate.cc)
target_link_libraries(chef_propagate  ${CHEF_LIBS} ${BOOST_LIBS})

add_python_extension(error_eater error_eater.cc)
target_link_libraries(error_eater ${BOOST_LIBS})

add_python_extension(mappers mappers.cc)
target_link_libraries(mappers  ${CHEF_LIBS} ${BOOST_LIBS})

add_python_extension(hist2d hist2d.cc)
target_link_libraries(hist2d ${BOOST_LIBS})
