link_libraries(${CHEF_LIBS} ${BOOST_LIBS})

#~ s2_containers.so s2_deposit.so s2_electric_field.so s2_solver_fftw.so macro_bunch_store.so 
#~ s2_diagnostics.so GaussSC.so populate.so

add_python_extension(GaussSC GaussSC.cc BasErs_field)

add_python_extension(s2_containers s2_containers_wrap.cc ltwt_containers.cc 
    mytimer.cc)
    
add_python_extension(s2_deposit s2_deposit_wrap.cc deposit.cc mytimer.cc)

add_python_extension(s2_solver_fftw s2_solver_fftw_wrap.cc solver_fftw.cc
    mytimer.cc fftw_helper.cc communicate.cc)
target_link_libraries(s2_solver_fftw ${FFTW2_R_MPI_LIB} ${FFTW2_MPI_LIB}
    ${FFTW2_R_LIB} ${FFTW2_LIB})

add_python_extension(s2_electric_field s2_electric_field_wrap.cc 
    electric_field.cc fftw_helper.cc mytimer.cc communicate.cc)
target_link_libraries(s2_electric_field ${FFTW2_R_MPI_LIB} ${FFTW2_MPI_LIB}
    ${FFTW2_R_LIB} ${FFTW2_LIB})
#~ target_link_libraries(s2_electric_field  s2_solver_fftw)
# s2_containers.so s2_solver_fftw.so

add_python_extension(s2_diagnostics s2_diagnostics_wrap.cc ltwt_containers.cc
    s2_diagnostics.cc)

add_python_extension(macro_bunch_store macro_bunch_store_wrap.cc 
    macro_bunch_store.cc ltwt_containers.cc)

add_python_extension(populate populate.cc populate_wrap.cc)
target_link_libraries(populate ${GSL_CONFIG_LIBS})