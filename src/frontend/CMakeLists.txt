configure_file("${SYNERGIA2_SOURCE_DIR}/src/frontend/synergia.in"
    "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" IMMEDIATE)

add_executable(synergia-resume synergia-resume.cc)
target_link_libraries(synergia-resume synergia_simulation synergia_command_line synergia_collective)

add_executable(synergia-archive-convert synergia-archive-convert.cc)
target_link_libraries(synergia-archive-convert synergia_simulation synergia_command_line synergia_collective)

install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" DESTINATION bin)
install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia-resume" DESTINATION bin)
install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia-archive-convert" DESTINATION bin)
