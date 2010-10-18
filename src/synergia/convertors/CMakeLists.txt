add_python_extension(pyconvertors convertors_wrap.cc)
target_link_libraries(pyconvertors parallel_utils ${Boost_LIBRARIES})

#add_subdirectory(tests)
