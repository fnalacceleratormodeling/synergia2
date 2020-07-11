add_executable(cxx_example cxx_example.cc)
target_link_libraries(cxx_example synergia_simulation synergia_collective synergia_serialization synergia_foundation synergia_parallel_utils synergia_bunch synergia_lattice)

copy_file(fodo.lat cxx_example)
copy_file (cxx_example.py cxx_example)

set(CLEAN_FILES
    cxx_covariance_matrix.xml
    cxx_lattice.xml
    cxx_means.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

