add_library(parallel_utils SHARED parallel_utils.cc commxx.cc)
add_python_extension(pyparallel_utils parallel_utils_wrap.cc)
target_link_libraries(pyparallel_utils parallel_utils ${Boost_LIBRARIES})

add_subdirectory(tests)
