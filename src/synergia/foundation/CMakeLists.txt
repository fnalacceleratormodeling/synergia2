file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(synergia_foundation four_momentum.cc reference_particle.cc
    distribution.cc multi_diagnostics.cc
    diagnostics_write_helper.cc)
target_link_libraries(synergia_foundation ${GSL_CONFIG_LIBS}
    synergia_parallel_utils synergia_hdf5_utils synergia_serialization
    ${Boost_SERIALIZATION_LIBRARIES})
if (BUILD_PYTHON_BINDINGS)
    add_python_extension(foundation foundation_wrap.cc)
    target_link_libraries(foundation synergia_foundation ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_foundation DESTINATION lib)
install(FILES
    distribution.h
    diagnostics_write_helper.h
    four_momentum.h
    math_constants.h
    multi_diagnostics.h
    physical_constants.h
    reference_particle.h
    DESTINATION include/synergia/foundation)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/foundation)
    install(TARGETS
        foundation
        DESTINATION lib/synergia/foundation)
endif ()

add_subdirectory(tests)
