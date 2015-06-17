add_python_extension(test_helper test_helper.cc)
target_link_libraries(test_helper ${Boost_LIBRARIES})

add_test(NAME test_numpy_multi_array_converter_py COMMAND 
    ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" 
        test_numpy_multi_array_converter.py)
