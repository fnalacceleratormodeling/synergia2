add_executable(gdfs1d gdfs1d.cc)
target_link_libraries(gdfs1d -lfftw3)

add_executable(fft_test fft_test.cc)
target_link_libraries(fft_test -lfftw3)

add_executable(test_slice test_slice.cc)
target_link_libraries(test_slice array_nd)

link_libraries(${CHEF_LIBS} -lfftw3 ${BOOST_LIBS} array_nd)

add_python_extension(GaussSC GaussSC.cc BasErs_field.cc)

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

add_python_extension(s2_diagnostics s2_diagnostics_wrap.cc ltwt_containers.cc
    s2_diagnostics.cc)

add_python_extension(macro_bunch_store macro_bunch_store_wrap.cc 
    macro_bunch_store.cc ltwt_containers.cc)

add_python_extension(populate populate.cc populate_wrap.cc)
target_link_libraries(populate ${GSL_CONFIG_LIBS})

add_python_extension(constraints constraints.cc constraints_wrap.cc)

add_python_extension(s2_solver_cylindrical cylindrical.cc cylindrical_wrap.cc
    field_domain.cc populate.cc ) # jfa fix this!!!!!!!!!!!!!!!!!!!!!!
target_link_libraries(s2_solver_cylindrical ${GSL_CONFIG_LIBS})
## jfa nb.: populate and gsl_config_libs are only part of s2_solver so that I can use the freakin' write to octave function for debugging. Much refoactoring necessary.
################
add_test(deposit ${SYNERGIA2_BINARY_DIR}/synergia-local test_deposit.py)
add_test(macro_bunch ${SYNERGIA2_BINARY_DIR}/synergia-local test_macro_bunch.py)
add_test(nd_array ${SYNERGIA2_BINARY_DIR}/synergia-local test_nd_array.py)
add_test(scalar_field ${SYNERGIA2_BINARY_DIR}/synergia-local test_scalar_field.py)
add_test(solver_fftw ${SYNERGIA2_BINARY_DIR}/synergia-local test_solver_fftw.py)
add_test(solver_cylindrical ${SYNERGIA2_BINARY_DIR}/synergia-local ${SYNERGIA2_SOURCE_DIR}/s2_fish/test_solver_cylindrical.py)
add_test(populate ${SYNERGIA2_BINARY_DIR}/synergia-local test_populate.py)
