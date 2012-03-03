configure_file("${SYNERGIA2_SOURCE_DIR}/src/frontend/synergia.in"
    "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" IMMEDIATE)

add_executable(synergia-resume synergia-resume.cc)
target_link_libraries(synergia-resume synergia_simulation synergia_command_line)

install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" DESTINATION bin)
install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia-resume" DESTINATION bin)
