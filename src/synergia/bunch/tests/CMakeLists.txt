add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_bunch test_bunch.cc)
target_link_libraries(test_bunch synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_bunch ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_bunch)
add_test(test_bunch_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_bunch.py)

add_test_executable(test_bunch_train test_bunch_train.cc)
target_link_libraries(test_bunch_train synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_bunch_train ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_bunch_train)
add_test(test_bunch_train_py
    ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()"
    test_bunch_train.py)

add_test_executable(test_bunch_train_mpi test_bunch_train_mpi.cc)
target_link_libraries(test_bunch_train_mpi synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_bunch_train_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_bunch_train_mpi)
add_test(test_bunch_train_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2  test_bunch_train_mpi)
add_test(test_bunch_train_mpi3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3  test_bunch_train_mpi)
add_test(test_bunch_train_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4  test_bunch_train_mpi)

#add_test_executable(test_bunch_with_diagnostics test_bunch_with_diagnostics.cc)
#target_link_libraries(test_bunch_with_diagnostics synergia_bunch
#    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#add_test(test_bunch_with_diagnostics ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_bunch_train)


add_custom_command(OUTPUT test_diagnostics_get_corr.icc
        test_diagnostics_get_emitx.icc
        test_diagnostics_get_emitxy.icc
        test_diagnostics_get_emitxyz.icc
        test_diagnostics_get_emity.icc
        test_diagnostics_get_emitz.icc
        test_diagnostics_get_mean.icc
        test_diagnostics_get_mom2.icc
        test_diagnostics_get_random_particles.icc
        test_diagnostics_get_std.icc
    COMMAND python test_diagnostics_crosscheck.py
    DEPENDS test_diagnostics_crosscheck.py)
add_test_executable(test_diagnostics test_diagnostics.cc test_diagnostics_get_corr.icc
    test_diagnostics_get_emitx.icc
    test_diagnostics_get_emitxy.icc
    test_diagnostics_get_emitxyz.icc
    test_diagnostics_get_emity.icc
    test_diagnostics_get_emitz.icc
    test_diagnostics_get_mean.icc
    test_diagnostics_get_mom2.icc
    test_diagnostics_get_random_particles.icc
    test_diagnostics_get_std.icc)
target_link_libraries(test_diagnostics synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics)
add_test(test_diagnostics_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_diagnostics.py)

add_test_executable(test_diagnostics_track test_diagnostics_track.cc)
target_link_libraries(test_diagnostics_track synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_track ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_track)
add_test(test_diagnostics_track_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_diagnostics_track.py)

add_test_executable(test_diagnostics_particles test_diagnostics_particles.cc)
add_dependencies(test_diagnostics_particles test_diagnostics_particles_get_mean.icc)
target_link_libraries(test_diagnostics_particles synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_particles ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_particles)
add_test(test_diagnostics_particles_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_diagnostics_particles.py)

add_test_executable(test_diagnostics_basic_mpi test_diagnostics_basic_mpi.cc)
target_link_libraries(test_diagnostics_basic_mpi synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_basic_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_basic_mpi)
add_test(test_diagnostics_basic_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_diagnostics_basic_mpi)
add_test(test_diagnostics_basic_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_diagnostics_basic_mpi)

add_test_executable(test_diagnostics_full2_mpi test_diagnostics_full2_mpi.cc)
target_link_libraries(test_diagnostics_full2_mpi synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_full2_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_full2_mpi)
add_test(test_diagnostics_full2_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_diagnostics_full2_mpi)
add_test(test_diagnostics_full2_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_diagnostics_full2_mpi)

add_test_executable(test_diagnostics_particles_mpi test_diagnostics_particles_mpi.cc)
target_link_libraries(test_diagnostics_particles_mpi synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_particles_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_particles_mpi)
add_test(test_diagnostics_particles_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_diagnostics_particles_mpi)
add_test(test_diagnostics_particles_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_diagnostics_particles_mpi)

add_test_executable(test_diagnostics_track_mpi test_diagnostics_track_mpi.cc)
target_link_libraries(test_diagnostics_track_mpi synergia_bunch
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_track_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_diagnostics_track_mpi)
add_test(test_diagnostics_track_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_diagnostics_track_mpi)
add_test(test_diagnostics_track_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_diagnostics_track_mpi)

#add_test_executable(test_multi_diagnostics_writer test_multi_diagnostics_write_helper.cc)
#target_link_libraries(test_multi_diagnostics_writer synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#add_test(test_multi_diagnostics_writer test_multi_diagnostics_writer)
#add_test(test_multi_diagnostics_writer_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_multi_diagnostics_writer.py)

add_test_executable(test_populate test_populate.cc)
target_link_libraries(test_populate synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_populate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_populate)
add_test(test_populate_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_populate.py)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
    "test_init_writers_full2.h5;test_init_writers.h5;test_init_writers_particles.h5;test_py_construct_full2.h5;test_py_construct.h5;test_py_get_diagnostics.h5;test_py_update_and_write_full2.h5;test_py_update_and_write.h5;test_py_write.h5;test_write_hdf5_full2.h5;test_write_hdf5.h5;test_write_hdf5_no_init.h5;test_write_hdf5_no_init_particles.h5;test_write_hdf5_particles.h5;test_writer_construct_full2.h5;test_writer_construct.h5;test_writer_get_count.h5;test_writer_get_diagnostics_sptr.h5;test_writer_is_dummy_false.h5;test_writer_update_and_write_full2.h5;test_writer_update_and_write.h5;test_writer_update_and_write_particles_0000.h5;test_writer_update_and_write_particles_0001.h5;test_writer_update_and_write_particles_0002.h5;test_writer_write.h5;dummy.h5;dummy1.h5;dummy2.h5")
