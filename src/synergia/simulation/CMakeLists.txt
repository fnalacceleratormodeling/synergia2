add_library(simulation SHARED operator.cc step.cc stepper.cc propagator.cc
    fast_mapping.cc chef_propagator.cc independent_operation.cc
    lattice_simulator.cc operation_extractor.cc)
target_link_libraries(simulation foundation bunch lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

add_python_extension(pysimulation simulation_wrap.cc)
target_link_libraries(pysimulation foundation bunch lattice simulation ${CHEF_LIBS} ${Boost_LIBRARIES})

install(TARGETS simulation DESTINATION lib)
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

add_subdirectory(tests)
