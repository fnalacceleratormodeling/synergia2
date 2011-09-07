add_library(synergia_simulation operator.cc step.cc stepper.cc
    propagate_actions.cc standard_diagnostics_actions.cc propagator.cc
    fast_mapping.cc chef_propagator.cc independent_operation.cc
    lattice_simulator.cc operation_extractor.cc aperture.cc)
target_link_libraries(synergia_simulation synergia_foundation synergia_bunch
    synergia_lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

if (BUILD_PYTHON_BINDINGS)
    add_python_extension(simulation simulation_wrap.cc)
    target_link_libraries(simulation synergia_foundation synergia_bunch
        synergia_lattice synergia_simulation
        ${CHEF_LIBS} ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_simulation DESTINATION lib)
install(FILES
    chef_propagator.h
    fast_mapping.h
    independent_operation.h
    lattice_simulator.h
    operation_extractor.h
    operator.h
    propagate_actions.h
    propagator.h
    step.h
    stepper.h
    DESTINATION include/synergia/simulation)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/simulation)
    install(TARGETS
        simulation
        DESTINATION lib/synergia/simulation)
endif ()

add_subdirectory(tests)
