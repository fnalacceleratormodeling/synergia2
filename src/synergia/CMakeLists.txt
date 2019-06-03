add_subdirectory(utils)
if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(convertors)
endif ()

add_subdirectory(foundation)
add_subdirectory(bunch)
add_subdirectory(lattice)
add_subdirectory(libFF)
add_subdirectory(simulation)
#add_subdirectory(optics)
#add_subdirectory(collective)
add_subdirectory(tools)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/synergia.cc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/synergia.cc" IMMEDIATE)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.py" IMMEDIATE)

add_library(synergia ${CMAKE_CURRENT_BINARY_DIR}/synergia.cc)
target_link_libraries(synergia synergia_parallel_utils synergia_hdf5_utils
    synergia_foundation synergia_bunch synergia_lattice synergia_simulation
    synergia_distributed_fft synergia_collective)

install(TARGETS synergia DESTINATION ${LIB_INSTALL_DIR})
install(FILES synergia.h DESTINATION ${INCLUDE_INSTALL_DIR}/synergia)
install(FILES
    __init__.py
    "${CMAKE_CURRENT_BINARY_DIR}/version.py"
    DESTINATION ${PYTHON_INSTALL_DIR}/synergia)
