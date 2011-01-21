add_library(synergia_simulation SHARED operator.cc step.cc stepper.cc propagator.cc
    fast_mapping.cc chef_propagator.cc independent_operation.cc
    lattice_simulator.cc operation_extractor.cc aperture.cc)
target_link_libraries(synergia_simulation synergia_foundation synergia_bunch
    synergia_lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

add_python_extension(simulation simulation_wrap.cc)
target_link_libraries(simulation synergia_foundation synergia_bunch
    synergia_lattice synergia_simulation
    ${CHEF_LIBS} ${Boost_LIBRARIES})

install(TARGETS synergia_simulation DESTINATION lib)
install(FILES
    chef_propagator.h
    fast_mapping.h
    independent_operation.h
    lattice_simulator.h
    operation_extractor.h
    operator.h
    propagator.h
    step.h
    stepper.h
    DESTINATION include/synergia/simulation)
install(FILES
    __init__.py
    DESTINATION lib/synergia/simulation)
install(TARGETS
    simulation
    DESTINATION lib/synergia/simulation)

add_subdirectory(tests)
