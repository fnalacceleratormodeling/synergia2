add_library(foundation SHARED four_momentum.cc reference_particle.cc)
#target_link_libraries(demo ${CHEF_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
add_subdirectory(tests)
