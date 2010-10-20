configure_file("${SYNERGIA2_SOURCE_DIR}/src/frontend/synergia.in"
    "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" IMMEDIATE)

install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synergia" DESTINATION bin)
