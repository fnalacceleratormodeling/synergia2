SET(COVERAGE_COMPILE_FLAGS "--coverage")
SET(COVERAGE_LINK_FLAGS    "--coverage")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")

add_executable(fodo_cxx fodo_cxx.cc)
target_link_libraries(fodo_cxx
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    )
