add_library(bunch SHARED bunch.cc)
#target_link_libraries(demo ${CHEF_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
add_subdirectory(tests)
