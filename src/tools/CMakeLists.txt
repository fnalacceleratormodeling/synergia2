configure_file("${SYNERGIA2_SOURCE_DIR}/src/tools/synmad8toxml.in"
    "${SYNERGIA2_BINARY_DIR}/src/frontend/synmad8toxml" IMMEDIATE)

install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/synmad8toxml" DESTINATION bin)

install(FILES
    mad8_to_xml.py
    DESTINATION lib/synergia_tools)
