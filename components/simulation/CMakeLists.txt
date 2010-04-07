add_library(simulation SHARED operator.cc step.cc stepper.cc propagator.cc
    fast_mapping.cc chef_propagator.cc independent_operation.cc
    independent_params.cc operation_extractor.cc)
target_link_libraries(simulation foundation bunch lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

#add_python_extension(pysimulation simulation_wrap.cc)
#target_link_libraries(pysimulation foundation lattice simulation ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

add_subdirectory(tests)
