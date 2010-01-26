add_executable(test_slice test_slice.cc)
target_link_libraries(test_slice array_nd)

# jfa: this is overkill
add_definitions(-DBOOST_TEST_DYN_LINK)

add_python_extension(GaussSC GaussSC.cc BasErs_field.cc)
target_link_libraries(GaussSC ${CHEF_LIBS} ${Boost_LIBRARIES})

add_python_extension(s2_containers s2_containers_wrap.cc ltwt_containers.cc 
    mytimer.cc)
target_link_libraries(s2_containers ${Boost_LIBRARIES})
    
add_python_extension(s2_deposit s2_deposit_wrap.cc deposit.cc mytimer.cc)
target_link_libraries(s2_deposit ${Boost_LIBRARIES})

add_python_extension(s2_solver_fftw s2_solver_fftw_wrap.cc solver_fftw.cc
    mytimer.cc fftw_helper.cc communicate.cc)
target_link_libraries(s2_solver_fftw ${Boost_LIBRARIES} ${PARALLEL_FFTW_LIBS})

add_python_extension(s2_electric_field s2_electric_field_wrap.cc 
    electric_field.cc fftw_helper.cc mytimer.cc communicate.cc)
target_link_libraries(s2_electric_field ${Boost_LIBRARIES} ${PARALLEL_FFTW_LIBS})

add_python_extension(s2_impedance_kick  impedance_kick.cc s2_impedance_kick_wrap.cc macro_bunch_store.cc)
target_link_libraries(s2_impedance_kick ${Boost_LIBRARIES})

add_python_extension(s2_diagnostics s2_diagnostics_wrap.cc ltwt_containers.cc
    s2_diagnostics.cc)
target_link_libraries(s2_diagnostics ${Boost_LIBRARIES})

add_python_extension(macro_bunch_store macro_bunch_store_wrap.cc 
    macro_bunch_store.cc ltwt_containers.cc)
target_link_libraries(macro_bunch_store  ${Boost_LIBRARIES})

add_python_extension(populate populate.cc populate_wrap.cc)
target_link_libraries(populate array_nd ${Boost_LIBRARIES} ${GSL_CONFIG_LIBS})

add_python_extension(constraints constraints.cc constraints_wrap.cc)
target_link_libraries(constraints array_nd ${Boost_LIBRARIES})

add_python_extension(s2_solver_transverse_fftw  solver_transverse_fftw.cc s2_solver_transverse_wrap.cc macro_bunch_store.cc)
target_link_libraries(s2_solver_transverse_fftw ${FFTW3_LIB} ${Boost_LIBRARIES})

add_python_extension(s2_solver_pplanes s2_solver_pplanes_wrap.cc  parallel_planes.cc parallelization.cc mytimer.cc)
target_link_libraries(s2_solver_pplanes ${Boost_LIBRARIES} ${FFTW3_LIB} ${PARALLEL_FFTW_LIBS}) 

add_python_extension(s2_solver_rectangular s2_solver_rectangular_wrap.cc  rectangular.cc  mytimer.cc)  
target_link_libraries(s2_solver_rectangular ${Boost_LIBRARIES})

add_library(s2_solver_cyl SHARED cylindrical.cc
    field_domain.cc populate.cc mytimer.cc) # jfa fix this!!!!!!!!!!!!!!!!!!!!!!
add_python_extension(s2_solver_cylindrical cylindrical_wrap.cc) # jfa fix this!!!!!!!!!!!!!!!!!!!!!!
target_link_libraries(s2_solver_cylindrical s2_solver_cyl array_nd ${FFTW3_LIB} ${Boost_LIBRARIES} ${GSL_CONFIG_LIBS})
## jfa nb.: populate and gsl_config_libs are only part of s2_solver so that I can use the freakin' write to octave function for debugging. Much refoactoring necessary.
################
add_test(deposit ${SYNERGIA2_BINARY_DIR}/synergia-local test_deposit.py)
add_test(macro_bunch ${SYNERGIA2_BINARY_DIR}/synergia-local test_macro_bunch.py)
add_test(nd_array ${SYNERGIA2_BINARY_DIR}/synergia-local test_nd_array.py)
add_test(scalar_field ${SYNERGIA2_BINARY_DIR}/synergia-local test_scalar_field.py)
add_test(solver_fftw ${SYNERGIA2_BINARY_DIR}/synergia-local test_solver_fftw.py)
add_test(solver_cylindrical ${SYNERGIA2_BINARY_DIR}/synergia-local ${SYNERGIA2_SOURCE_DIR}/s2_fish/test_solver_cylindrical.py)
add_test(populate ${SYNERGIA2_BINARY_DIR}/synergia-local test_populate.py)
