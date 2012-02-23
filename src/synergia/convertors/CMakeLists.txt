add_python_extension(convertors convertors_wrap.cc)
target_link_libraries(convertors synergia_parallel_utils synergia_serialization
    ${Boost_LIBRARIES})

install(FILES
    __init__.py
    DESTINATION lib/synergia/convertors)
install(TARGETS
    convertors
    DESTINATION lib/synergia/convertors)

#add_subdirectory(tests)
