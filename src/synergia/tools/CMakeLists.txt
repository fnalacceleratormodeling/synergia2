file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/three_bump.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(FILES
    __init__.py
    three_bump.py
    DESTINATION lib/synergia/tools)

add_subdirectory(tests)
