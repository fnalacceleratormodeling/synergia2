set(tools
    synmad8toxml
    syninspecth5
    syndiagplot
    syntrackplot
    synpoincareplot
    synbeamplot
    synlatticeview
    synprintmaps
    synpart2txt
    syntrack2txt)

foreach(tool ${tools})
    configure_file("${SYNERGIA2_SOURCE_DIR}/src/tools/${tool}.in"
        "${SYNERGIA2_BINARY_DIR}/src/tools/${tool}" IMMEDIATE)
    install(PROGRAMS "${SYNERGIA2_BINARY_DIR}/src/tools/${tool}" DESTINATION bin)
endforeach(tool tools)

install(FILES
    mad8_to_xml.py
    inspect_h5.py
    diag_plot.py
    track_plot.py
    poincare_plot.py
    beam_plot.py
    lattice_view.py
    print_maps.py
    particles_to_text.py
    tracks_to_text.py
    DESTINATION lib/synergia_tools)
