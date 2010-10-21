add_library(parallel_utils SHARED parallel_utils.cc commxx.cc)
add_python_extension(pyparallel_utils parallel_utils_wrap.cc)
target_link_libraries(pyparallel_utils parallel_utils ${Boost_LIBRARIES})

add_library(hdf5_utils SHARED hdf5_writer.cc hdf5_serial_writer.cc
    hdf5_chunked_array2d_writer.cc)
target_link_libraries(hdf5_utils ${HDF5_LIBRARIES})

add_library(distributed_fft SHARED distributed_fft3d.cc)
target_link_libraries(distributed_fft ${PARALLEL_FFTW_LIBRARIES} parallel_utils)

install(TARGETS parallel_utils hdf5_utils DESTINATION lib)
install(FILES
    boost_test_mpi_fixture.h
    comm_converter.h
    commxx.h
    container_conversions.h
    distributed_fft3d.h
    floating_point.h
    hdf5_chunked_array2d_writer.h
    hdf5_misc.h
    hdf5_serial_writer.h
    hdf5_serial_writer.tcc
    hdf5_writer.h
    hdf5_writer.tcc
    multi_array_check_equal.h
    multi_array_print.h
    multi_array_serialization.h
    multi_array_typedefs.h
    numpy_multi_ref_converter.h
    parallel_utils.h
    xml_serialization.h
    DESTINATION include/synergia/utils)

add_subdirectory(tests)
