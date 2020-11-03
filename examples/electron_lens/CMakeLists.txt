add_executable(cxx_offdiag s3_cxx_offdiag.cc)
target_link_libraries(cxx_offdiag
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    )

copy_file(adjusted_lattice.lsx cxx_offdiag)

