add_subdirectory(utils)
if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(convertors)
endif ()
add_subdirectory(foundation)
add_subdirectory(bunch)
add_subdirectory(lattice)
add_subdirectory(simulation)
add_subdirectory(optics)
add_subdirectory(collective)
add_subdirectory(tools)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file("synergia.cc.in" "synergia.cc")

configure_file("version.py.in" "version.py")

add_library(synergia ${CMAKE_CURRENT_BINARY_DIR}/synergia.cc)
target_link_libraries(synergia synergia_parallel_utils synergia_hdf5_utils
    synergia_foundation synergia_bunch synergia_lattice synergia_simulation
    synergia_distributed_fft synergia_collective)

install(TARGETS synergia DESTINATION lib)
install(FILES synergia.h DESTINATION include/synergia)
install(FILES
    __init__.py
    "${CMAKE_CURRENT_BINARY_DIR}/version.py"
    DESTINATION lib/synergia)
