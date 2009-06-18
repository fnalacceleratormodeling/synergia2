configure_file("${SYNERGIA2_SOURCE_DIR}/frontend/synergia.in" 
    "${SYNERGIA2_BINARY_DIR}/frontend/synergia" IMMEDIATE)

install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/frontend/synergia" DESTINATION bin)
