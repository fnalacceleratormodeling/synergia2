link_libraries(${CHEF_LIBS} ${BOOST_LIBS})

add_python_extension(GaussSC GaussSC.cc BasErs_field)

add_python_extension(s2_containers s2_containers_wrap.cc ltwt_containers.cc 
    mytimer.cc)

add_python_extension(macro_bunch_store macro_bunch_store_wrap.cc 
    macro_bunch_store.cc ltwt_containers.cc)

add_python_extension(populate populate.cc)
