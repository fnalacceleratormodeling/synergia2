set(tools
    synmad8toxml
    syninspecth5
    syndiagplot
    synbeamplot)

foreach(tool ${tools})
    configure_file("${SYNERGIA2_SOURCE_DIR}/src/tools/${tool}.in"
        "${SYNERGIA2_BINARY_DIR}/src/frontend/${tool}" IMMEDIATE)
    install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/frontend/${tool}" DESTINATION bin)
endforeach(tool tools)

install(FILES
    mad8_to_xml.py
    inspect_h5.py
    diag_plot.py
    beam_plot.py
    DESTINATION lib/synergia_tools)
