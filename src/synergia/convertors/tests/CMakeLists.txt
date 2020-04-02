add_python_extension(test_helper test_helper.cc)
target_link_libraries(test_helper ${Boost_LIBRARIES})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_numpy_multi_array_converter.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_numpy_multi_ref_converter.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_numpy_multi_array_converter_py COMMAND 
    ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" 
        test_numpy_multi_array_converter.py)

add_test(NAME test_numpy_multi_ref_converter_py COMMAND 
    ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" 
        test_numpy_multi_ref_converter.py)

