SET(CMAKE_CXX_COMPILER ${mpicxx})

add_library(GaussSC MODULE GaussSC.cc BasErs_field)
target_link_libraries(GaussSC  ${CHEF_LIBS} ${BOOST_LIBS})
set_target_properties(GaussSC PROPERTIES PREFIX "")

add_library(s2_containers MODULE s2_containers_wrap.cc ltwt_containers.cc mytimer.cc)
target_link_libraries(s2_containers ${CHEF_LIBS} ${BOOST_LIBS})
set_target_properties(s2_containers PROPERTIES PREFIX "")

