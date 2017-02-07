set(tools
    synmad8toxml
    synmad8tolsx
    syninspecth5
    syndiagplot
    syntrackplot
    synpoincareplot
    synbeamplot
    synlatticeview
    synlatticefns
    synprintmaps
    synpart2txt
    syntrack2txt
    synpart2madxtxt
    syndecodemad8twiss)

foreach(tool ${tools})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${tool}.in"
        "${CMAKE_CURRENT_BINARY_DIR}/${tool}" IMMEDIATE)
    install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${tool}"
        DESTINATION ${BIN_INSTALL_DIR})
endforeach(tool tools)

install(FILES
    mad8_to_xml.py
    mad8_to_lsx.py
    inspect_h5.py
    diag_plot.py
    track_plot.py
    poincare_plot.py
    beam_plot.py
    lattice_view.py
    lattice_fns.py
    print_maps.py
    particles_to_text.py
    tracks_to_text.py
    particles_to_madX_text.py
    decode_mad8_twiss.py
    DESTINATION ${PYTHON_INSTALL_DIR}/synergia_tools)
