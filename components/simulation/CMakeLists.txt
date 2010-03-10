add_library(simulation SHARED stepper.cc)
target_link_libraries(simulation foundation lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

#add_python_extension(pysimulation simulation_wrap.cc)
#target_link_libraries(pysimulation foundation lattice simulation ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

add_subdirectory(tests)
