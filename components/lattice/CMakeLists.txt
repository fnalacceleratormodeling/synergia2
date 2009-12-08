#include_directories(/usr/include/glib-2.0 /usr/lib64/glib-2.0/include)
#add_executable(demo demo.cc lattice_element.cc lattice.cc four_momentum.cc reference_particle.cc)
#target_link_libraries(demo ${CHEF_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

add_subdirectory(tests)
