add_executable(foborodobo32 foborodobo32.cc)
target_link_libraries(foborodobo32 basic_toolkit beamline mxyzptlk physics_toolkit synergia_foundation synergia_lattice synergia_parallel_utils synergia_simulation synergia_collective synergia_serialization synergia_hdf5_utils synergia_distributed_fft synergia)

#set(CLEAN_FILES
#    cxx_covariance_matrix.xml
#    cxx_lattice.xml
#    cxx_means.xml)

#set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

