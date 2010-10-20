add_executable(cxx_example cxx_example.cc)
target_link_libraries(cxx_example simulation)

set(CLEAN_FILES
    cxx_covariance_matrix.xml
    cxx_lattice.xml
    cxx_means.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

