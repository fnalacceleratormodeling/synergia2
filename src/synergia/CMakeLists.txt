add_subdirectory(utils)
add_subdirectory(convertors)
add_subdirectory(foundation)
add_subdirectory(bunch)
add_subdirectory(lattice)
add_subdirectory(simulation)
add_subdirectory(optics)
add_subdirectory(collective)

configure_file("${SYNERGIA2_SOURCE_DIR}/src/synergia/synergia.cc.in"
    "${SYNERGIA2_BINARY_DIR}/src/synergia/synergia.cc" IMMEDIATE)

configure_file("${SYNERGIA2_SOURCE_DIR}/src/synergia/version.py.in"
    "${SYNERGIA2_BINARY_DIR}/src/synergia/version.py" IMMEDIATE)

add_library(synergia SHARED synergia.cc)
target_link_libraries(synergia synergia_parallel_utils synergia_hdf5_utils
    synergia_foundation synergia_bunch synergia_lattice synergia_simulation
    synergia_distributed_fft synergia_collective)

install(TARGETS synergia DESTINATION lib)
install(FILES synergia.h DESTINATION include/synergia)
install(FILES
    __init__.py
    version.py
    DESTINATION lib/synergia)
