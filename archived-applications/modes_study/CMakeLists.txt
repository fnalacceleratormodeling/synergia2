add_executable(cxx_modes cxx_modes.cc excite_mode_actions.cc)
target_link_libraries(cxx_modes synergia_simulation synergia_collective synergia_bunch)

add_executable(cxx_modes_resume cxx_modes_resume.cc excite_mode_actions.cc)
target_link_libraries(cxx_modes_resume synergia_simulation synergia_collective synergia_bunch)



set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

