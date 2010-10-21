add_library(foundation SHARED four_momentum.cc reference_particle.cc distribution.cc)
target_link_libraries(foundation ${GSL_CONFIG_LIBS} parallel_utils ${Boost_SERIALIZATION_LIBRARIES})
add_python_extension(pyfoundation foundation_wrap.cc)
target_link_libraries(pyfoundation foundation ${Boost_LIBRARIES})

install(TARGETS foundation DESTINATION lib)
install(FILES
    distribution.h
    four_momentum.h
    math_constants.h
    physical_constants.h
    reference_particle.h
    DESTINATION include/synergia/foundation)

add_subdirectory(tests)
