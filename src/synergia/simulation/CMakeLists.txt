file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(synergia_simulation operator.cc step.cc stepper.cc
    independent_stepper.cc independent_stepper_elements.cc
    split_operator_stepper.cc split_operator_stepper_elements.cc
    split_operator_stepper_choice.cc
    propagate_actions.cc lattice_elements_actions.cc diagnostics_actions.cc propagator.cc
    fast_mapping.cc  fast_normal_form.cc chef_propagator.cc independent_operation.cc
    lattice_simulator.cc operation_extractor.cc aperture.cc
    aperture_operation.cc aperture_operation_extractor.cc populate_stationary.cc bunch_simulator.cc
    bunch_train_simulator.cc resume.cc dense_mapping.cc calculate_closed_orbit.cc)
target_link_libraries(synergia_simulation synergia_foundation synergia_bunch
    synergia_lattice ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

if (BUILD_PYTHON_BINDINGS)
    add_python_extension(simulation simulation_wrap.cc)
    target_link_libraries(simulation synergia_foundation synergia_bunch
        synergia_lattice synergia_simulation
        ${CHEF_LIBS} ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_simulation DESTINATION lib)
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
    DESTINATION include/synergia/simulation)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/simulation)
    install(TARGETS
        simulation
        DESTINATION lib/synergia/simulation)
endif ()

add_subdirectory(tests)
