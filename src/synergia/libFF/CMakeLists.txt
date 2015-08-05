add_library(synergia_libFF
    ff_element.cc
    ff_drift.cc
    ff_sbend.cc
    ff_rbend.cc
    ff_multipole.cc
    ff_quadrupole.cc
    ff_sextupole.cc
    ff_octupole.cc
    ff_rfcavity.cc
    ff_kicker.cc
    ff_hkicker.cc
    ff_vkicker.cc
    ff_element_map.cc)
target_link_libraries(synergia_libFF
   synergia_bunch synergia_lattice synergia_serialization)

#if (BUILD_PYTHON_BINDINGS)
#    add_python_extension(foundation foundation_wrap.cc)
#    target_link_libraries(foundation synergia_foundation ${Boost_LIBRARIES})
#endif ()

install(TARGETS synergia_libFF DESTINATION lib)
install(FILES
    ff_element.h
    DESTINATION include/synergia/libFF)

add_subdirectory(tests)
