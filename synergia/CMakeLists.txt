include_directories("..")

add_python_extension(apply_map apply_map.cc)
target_link_libraries(apply_map  ${CHEF_LIBS} ${Boost_LIBRARIES})

add_python_extension(chef_propagate chef_propagate.cc)
target_link_libraries(chef_propagate  ${CHEF_LIBS} ${Boost_LIBRARIES})

add_python_extension(error_eater error_eater.cc)
target_link_libraries(error_eater ${Boost_LIBRARIES})

add_python_extension(mappers mappers.cc)
target_link_libraries(mappers  ${CHEF_LIBS} ${Boost_LIBRARIES})

add_python_extension(hist2d hist2d.cc)
target_link_libraries(hist2d ${Boost_LIBRARIES})
