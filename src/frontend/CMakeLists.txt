configure_file("${CMAKE_CURRENT_SOURCE_DIR}/synergia.in"
    "${CMAKE_CURRENT_BINARY_DIR}/synergia" IMMEDIATE)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/synergia-pyresume.in"
    "${CMAKE_CURRENT_BINARY_DIR}/synergia-pyresume" IMMEDIATE)

add_executable(synergia-resume synergia-resume.cc)
target_link_libraries(synergia-resume synergia_simulation synergia_command_line synergia_collective)

add_executable(synergia-archive-convert synergia-archive-convert.cc)
target_link_libraries(synergia-archive-convert synergia_simulation synergia_command_line synergia_collective)

install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/synergia"
        DESTINATION ${BIN_INSTALL_DIR})
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/synergia-resume"
        DESTINATION ${BIN_INSTALL_DIR})
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/synergia-pyresume"
        DESTINATION ${BIN_INSTALL_DIR})
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/synergia-archive-convert"
        DESTINATION ${BIN_INSTALL_DIR})

install(FILES
    synergia_resume.py
    DESTINATION ${PYTHON_INSTALL_DIR}/synergia_tools)
