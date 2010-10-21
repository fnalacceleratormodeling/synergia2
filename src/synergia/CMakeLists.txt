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
target_link_libraries(synergia parallel_utils hdf5_utils foundation bunch
    lattice simulation distributed_fft space_charge)

install(TARGETS synergia DESTINATION lib)
install(FILES synergia.h DESTINATION include/synergia)
