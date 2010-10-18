add_subdirectory(utils)
add_subdirectory(convertors)
add_subdirectory(foundation)
add_subdirectory(bunch)
add_subdirectory(lattice)
add_subdirectory(simulation)
add_subdirectory(optics)
add_subdirectory(collective)

add_executable(cxx_example cxx_example.cc)
target_link_libraries(cxx_example simulation)

set(CLEAN_FILES
    cxx_covariance_matrix.xml
    cxx_lattice.xml
    cxx_means.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

