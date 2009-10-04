link_directories(..)
link_directories(../../../utils)

add_library(bunch SHARED bunch.cc)
target_link_libraries(bunch parallel_utils)
#target_link_libraries(demo ${CHEF_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
add_subdirectory(tests)
