add_executable(print_one_turn_map print_one_turn_map.cc)
add_executable(normal_form_example normal_form_example.cc)
add_executable(normal_form_mi normal_form_mi.cc)
target_link_libraries(print_one_turn_map synergia_simulation synergia_collective)
target_link_libraries(normal_form_example synergia_simulation synergia_collective)
target_link_libraries(normal_form_mi synergia_simulation synergia_collective)

set(CLEAN_FILES
    cxx_covariance_matrix.xml
    cxx_lattice.xml
    cxx_means.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

