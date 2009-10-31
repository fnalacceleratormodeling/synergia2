link_directories(.. ../../../utils ../foundation)

add_library(bunch SHARED fixed_t_z_converter.cc bunch.cc diagnostics.cc)
target_link_libraries(bunch parallel_utils)
add_python_extension(pybunch bunch_wrap.cc)
target_link_libraries(pybunch bunch parallel_utils foundation ${Boost_LIBRARIES}) 

add_subdirectory(doc)
add_subdirectory(tests)
