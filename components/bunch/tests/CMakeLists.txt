include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_bunch test_bunch.cc)
target_link_libraries(test_bunch bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_bunch test_bunch)
add_test(test_bunch_py nosetests test_bunch.py)

add_custom_command(OUTPUT test_diagnostics_get_mean.icc
    COMMAND python test_diagnostics_crosscheck.py
    DEPENDS test_diagnostics_crosscheck.py)
add_executable(test_diagnostics test_diagnostics.cc)
add_dependencies(test_diagnostics test_diagnostics_get_mean.icc)
target_link_libraries(test_diagnostics bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics test_diagnostics)
add_test(test_diagnostics_py nosetests test_diagnostics.py)

add_executable(test_diagnostics_writer test_diagnostics_writer.cc)
target_link_libraries(test_diagnostics_writer bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics_writer test_diagnostics_writer)
add_test(test_diagnostics_writer_py nosetests test_diagnostics_writer.py)

add_executable(test_populate test_populate.cc)
target_link_libraries(test_populate bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_populate test_populate)
add_test(test_populate_py nosetests test_populate.py)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
    "test_full2.h5;test_init_writers_full2.h5;test_init_writers.h5;test_py_construct_full2.h5;test_py_construct.h5;test_py_get_diagnostics.h5;test_py_update_and_write_full2.h5;test_py_update_and_write.h5;test_py_write.h5;test_write_hdf5_full2.h5;test_write_hdf5.h5;test_write_hdf5_no_init.h5;test_writer_construct_full2.h5;test_writer_construct.h5;test_writer_get_diagnostics_sptr.h5;test_writer.h5;test_writer_update_and_write_full2.h5;test_writer_update_and_write.h5;test_writer_write.h5;")
