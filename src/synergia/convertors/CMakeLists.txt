file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_python_extension(convertors convertors_wrap.cc)
target_link_libraries(convertors synergia_parallel_utils synergia_serialization
    ${Boost_LIBRARIES})

install(FILES
    __init__.py
    DESTINATION ${PYTHON_INSTALL_DIR}/synergia/convertors)
install(TARGETS
    convertors
    DESTINATION ${PYTHON_INSTALL_DIR}/synergia/convertors)

add_subdirectory(tests)
