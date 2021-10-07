if(${Kokkos_ENABLE_CUDA})
    set_source_files_properties(checkpoint_json.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(lattice_simulator_host.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()


add_library(synergia_simulation 
    operator.cc 
    step.cc 
    stepper.cc
    #independent_stepper.cc 
    independent_stepper_elements.cc
    split_operator_stepper.cc 
    split_operator_stepper_elements.cc
    #split_operator_stepper_choice.cc
    #propagate_actions.cc 
    #lattice_elements_actions.cc 
    #diagnostics_actions.cc 
    propagator.cc
    #fast_mapping.cc  
    #fast_normal_form.cc 
    #chef_propagator.cc 
    independent_operation.cc
    lattice_simulator.cc 
    lattice_simulator_host.cc 
    operation_extractor.cc 
    #aperture.cc
    aperture_operation.cc 
    #aperture_operation_extractor.cc 
    #populate_stationary.cc 
    bunch_simulator.cc
    #bunch_train_simulator.cc 
    #resume.cc 
    #dense_mapping.cc 
    #calculate_closed_orbit.cc 
    #diagnostics_normal_form.cc
    checkpoint.cc
    checkpoint_json.cc
    )

target_link_libraries(synergia_simulation 
    synergia_foundation 
    synergia_bunch
    synergia_lattice 
    synergia_collective
    )

target_link_options(synergia_simulation
    PRIVATE ${LINKER_OPTIONS}
    )

if (BUILD_PYTHON_BINDINGS)

    pybind11_add_module(simulation MODULE NO_EXTRAS simulation_pywrap.cc)
    target_link_libraries(simulation PRIVATE synergia_simulation synergia_pydiag)

    copy_file(__init__.py simulation)

endif ()

install(TARGETS synergia_simulation 
    DESTINATION ${LIB_INSTALL_DIR})

install(FILES
    bunch_simulator.h
    bunch_train_simulator.h
    chef_propagator.h
    diagnostics_actions.h
    fast_mapping.h
    fast_mapping.tcc
    fast_normal_form.h
    independent_operation.h
    aperture_operation.h
    aperture_operation.tcc
    aperture_operation_extractor.h
    lattice_simulator.h
    operation_extractor.h
    operator.h
    propagate_actions.h
    lattice_elements_actions.h
    propagator.h
    resume.h
    step.h
    stepper.h
    independent_stepper.h
    independent_stepper_elements.h
    split_operator_stepper.h
    split_operator_stepper_elements.h
    split_operator_stepper_choice.h
    populate_stationary.h
    dense_mapping.h
    calculate_closed_orbit.h
    diagnostics_normal_form.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/simulation
    )

#if (BUILD_PYTHON_BINDINGS)
#    install(FILES
#        __init__.py
#        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/simulation)
#    install(TARGETS
#        simulation
#        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/simulation)
#endif ()

add_subdirectory(tests)
