add_library(foundation SHARED four_momentum.cc reference_particle.cc distribution.cc)
target_link_libraries(foundation ${GSL_CONFIG_LIBS} parallel_utils)
add_python_extension(pyfoundation foundation_wrap.cc)
target_link_libraries(pyfoundation foundation ${Boost_LIBRARIES})

add_subdirectory(doc)
add_subdirectory(tests)
