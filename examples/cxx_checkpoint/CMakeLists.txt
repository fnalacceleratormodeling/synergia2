add_executable(cxx_checkpoint cxx_checkpoint.cc)
target_link_libraries(cxx_checkpoint synergia_simulation synergia_collective)

set(CLEAN_FILES
    cxx_covariance_matrix.xml
    cxx_lattice.xml
    cxx_means.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

