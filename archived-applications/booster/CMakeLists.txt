add_executable(cxx_booster cxx_booster.cc)
target_link_libraries(cxx_booster synergia_simulation synergia_collective)

add_executable(cxx_booster_resume cxx_booster_resume.cc)
target_link_libraries(cxx_booster_resume synergia_simulation synergia_collective)

set(CLEAN_FILES
    booster_lattice.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
