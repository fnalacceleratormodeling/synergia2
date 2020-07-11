file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(synergia_foundation four_momentum.cc reference_particle.cc
    distribution.cc multi_diagnostics.cc
    diagnostics_write_helper.cc)
#target_include_directories()
target_link_libraries(synergia_foundation 
    synergia_parallel_utils synergia_hdf5_utils synergia_serialization
    lsexpr
    ${GSL_LIBRARY}
    ${GSL_CBLAS_LIBRARY}
    ${Boost_SERIALIZATION_LIBRARIES})

add_python_extension(foundation foundation_wrap.cc)
target_link_libraries(foundation synergia_foundation ${Boost_LIBRARIES})

install(TARGETS synergia_foundation DESTINATION ${LIB_INSTALL_DIR})
install(FILES
    distribution.h
    diagnostics_write_helper.h
    four_momentum.h
    math_constants.h
    multi_diagnostics.h
    physical_constants.h
    reference_particle.h
    DESTINATION include/synergia/foundation)

install(FILES
        __init__.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)
install(TARGETS
        foundation
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)

add_subdirectory(tests)
