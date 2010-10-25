add_python_extension(convertors convertors_wrap.cc)
target_link_libraries(convertors synergia_parallel_utils ${Boost_LIBRARIES})

#add_subdirectory(tests)
