add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_rectangular_grid_domain test_rectangular_grid_domain.cc)
target_link_libraries(test_rectangular_grid_domain synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_rectangular_grid_domain test_rectangular_grid_domain)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_deposit test_deposit.cc)
target_link_libraries(test_deposit synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_deposit COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_deposit>)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_deposit_xyz test_deposit_xyz.cc)
target_link_libraries(test_deposit_xyz synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_deposit_xyz COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_deposit_xyz>)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)


add_test_executable(test_rectangular_grid test_rectangular_grid.cc)
target_link_libraries(test_rectangular_grid synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_rectangular_grid test_rectangular_grid)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_distributed_rectangular_grid test_distributed_rectangular_grid.cc)
target_link_libraries(test_distributed_rectangular_grid synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_rectangular_grid test_distributed_rectangular_grid)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_distributed_rectangular_grid_mpi test_distributed_rectangular_grid_mpi.cc)
target_link_libraries(test_distributed_rectangular_grid_mpi synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_distributed_rectangular_grid_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_distributed_rectangular_grid_mpi>)
add_test(NAME test_distributed_rectangular_grid_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_distributed_rectangular_grid_mpi>)
add_test(NAME test_distributed_rectangular_grid_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_distributed_rectangular_grid_mpi>)
add_test(NAME test_distributed_rectangular_grid_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_distributed_rectangular_grid_mpi>)
#add_test(NAME test_bunch_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_space_charge_rectangular
    test_space_charge_rectangular.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_rectangular synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_rectangular COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_rectangular>)
add_test(NAME test_space_charge_rectangular_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_space_charge_rectangular>)
add_test(NAME test_space_charge_rectangular_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_space_charge_rectangular>)
add_test(NAME test_space_charge_rectangular_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_space_charge_rectangular>)

add_test_executable(test_interpolate_rectangular_zyx
    test_interpolate_rectangular_zyx.cc gaussian_charge_density.cc)
target_link_libraries(test_interpolate_rectangular_zyx synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_interpolate_rectangular_zyx test_interpolate_rectangular_zyx)

add_test_executable(test_interpolate_rectangular_xyz
    test_interpolate_rectangular_xyz.cc gaussian_charge_density.cc)
target_link_libraries(test_interpolate_rectangular_xyz synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_interpolate_rectangular_xyz test_interpolate_rectangular_xyz)


add_test_executable(test_space_charge_3d_open_hockney
    test_space_charge_3d_open_hockney.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney>)

add_test_executable(test_space_charge_3d_open_hockney_mpi
    test_space_charge_3d_open_hockney_mpi.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney_mpi synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney_mpi>)
add_test(NAME test_space_charge_3d_open_hockney_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_space_charge_3d_open_hockney_mpi>)
add_test(NAME test_space_charge_3d_open_hockney_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_space_charge_3d_open_hockney_mpi>)
add_test(NAME test_space_charge_3d_open_hockney_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_space_charge_3d_open_hockney_mpi>)

add_test_executable(test_space_charge_3d_open_hockney2_mpi
    test_space_charge_3d_open_hockney2_mpi.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney2_mpi synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney2_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney2_mpi>)
add_test(NAME test_space_charge_3d_open_hockney2_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_space_charge_3d_open_hockney2_mpi>)
add_test(NAME test_space_charge_3d_open_hockney2_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_space_charge_3d_open_hockney2_mpi>)
add_test(NAME test_space_charge_3d_open_hockney2_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_space_charge_3d_open_hockney2_mpi>)

add_test_executable(test_space_charge_3d_open_hockney2
    test_space_charge_3d_open_hockney2.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney2 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney2>)

add_test_executable(test_space_charge_3d_open_hockney3
    test_space_charge_3d_open_hockney3.cc linear_cylindrical_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney3 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney3>)

add_test_executable(test_space_charge_3d_open_hockney4
    test_space_charge_3d_open_hockney4.cc uniform_cylindrical_charge_density.cc)
target_link_libraries(test_space_charge_3d_open_hockney4 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_3d_open_hockney4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_3d_open_hockney4>)

add_test_executable(test_space_charge_2d_open_hockney
    test_space_charge_2d_open_hockney.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney>)

add_test_executable(test_space_charge_2d_open_hockney_mpi
    test_space_charge_2d_open_hockney_mpi.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney_mpi
    synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney_mpi>)
add_test(NAME test_space_charge_2d_open_hockney_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_space_charge_2d_open_hockney_mpi>)
add_test(NAME test_space_charge_2d_open_hockney_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_space_charge_2d_open_hockney_mpi>)
add_test(NAME test_space_charge_2d_open_hockney_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_space_charge_2d_open_hockney_mpi>)

add_test_executable(test_space_charge_2d_open_hockney2
    test_space_charge_2d_open_hockney2.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney2 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney2>)

add_test_executable(test_space_charge_2d_open_hockney2_mpi
    test_space_charge_2d_open_hockney2_mpi.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney2_mpi
    synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney2_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney2_mpi>)
add_test(NAME test_space_charge_2d_open_hockney2_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_space_charge_2d_open_hockney2_mpi>)
add_test(NAME test_space_charge_2d_open_hockney2_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_space_charge_2d_open_hockney2_mpi>)
add_test(NAME test_space_charge_2d_open_hockney2_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_space_charge_2d_open_hockney2_mpi>)

add_test_executable(test_space_charge_2d_open_hockney3
    test_space_charge_2d_open_hockney3.cc linear_cylindrical_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney3 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney3>)

add_test_executable(test_space_charge_2d_open_hockney4
    test_space_charge_2d_open_hockney4.cc uniform_cylindrical_charge_density.cc)
target_link_libraries(test_space_charge_2d_open_hockney4 synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_open_hockney4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_open_hockney4>)

add_test_executable(test_space_charge_2d_bassetti_erskine
test_space_charge_2d_bassetti_erskine.cc gaussian_charge_density.cc)
target_link_libraries(test_space_charge_2d_bassetti_erskine synergia_collective
    synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_space_charge_2d_bassetti_erskine COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_space_charge_2d_bassetti_erskine>)

add_test_executable(test_impedance test_impedance.cc)
target_link_libraries(test_impedance synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_impedance COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_impedance>)
add_test(NAME test_impedance_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_impedance>)
add_test(NAME test_impedance_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_impedance>)
add_test(NAME test_impedance_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_impedance>)

add_test_executable(test_wake_field test_wake_field.cc)
target_link_libraries(test_wake_field synergia_collective
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_wake_field COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_wake_field>)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/wake_files
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
#    "test_init_writers_full2.h5;test_init_writers.h5;test_init_writers_particles.h5;test_py_construct_full2.h5;test_py_construct.h5;test_py_get_diagnostics.h5;test_py_update_and_write_full2.h5;test_py_update_and_write.h5;test_py_write.h5;test_write_hdf5_full2.h5;test_write_hdf5.h5;test_write_hdf5_no_init.h5;test_write_hdf5_no_init_particles.h5;test_write_hdf5_particles.h5;test_writer_construct_full2.h5;test_writer_construct.h5;test_writer_get_count.h5;test_writer_get_diagnostics_sptr.h5;test_writer_is_dummy_false.h5;test_writer_update_and_write_full2.h5;test_writer_update_and_write.h5;test_writer_update_and_write_particles_0000.h5;test_writer_update_and_write_particles_0001.h5;test_writer_update_and_write_particles_0002.h5;test_writer_write.h5")
