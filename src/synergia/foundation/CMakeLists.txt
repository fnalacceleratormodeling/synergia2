add_library(synergia_foundation SHARED four_momentum.cc reference_particle.cc distribution.cc)
target_link_libraries(synergia_foundation ${GSL_CONFIG_LIBS} synergia_parallel_utils ${Boost_SERIALIZATION_LIBRARIES})
add_python_extension(foundation foundation_wrap.cc)
target_link_libraries(foundation synergia_foundation ${Boost_LIBRARIES})

install(TARGETS synergia_foundation DESTINATION lib)
install(FILES
    distribution.h
    four_momentum.h
    math_constants.h
    physical_constants.h
    reference_particle.h
    DESTINATION include/synergia/foundation)

add_subdirectory(tests)
