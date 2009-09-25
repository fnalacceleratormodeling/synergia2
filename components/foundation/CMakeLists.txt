add_library(foundation SHARED four_momentum.cc reference_particle.cc)
add_python_extension(pyfoundation foundation_wrap.cc)
target_link_libraries(pyfoundation foundation ${Boost_LIBRARIES}) 

add_subdirectory(doc)
add_subdirectory(tests)
