add_library(parallel_utils SHARED parallel_utils.cc commxx.cc)
add_python_extension(pyparallel_utils parallel_utils_wrap.cc)
target_link_libraries(pyparallel_utils parallel_utils ${Boost_LIBRARIES})

add_library(hdf5_utils SHARED hdf5_writer.cc hdf5_serial_writer.cc
    hdf5_chunked_array2d_writer.cc)
target_link_libraries(hdf5_utils ${HDF5_LIBRARIES})

add_library(distrbuted_fft SHARED distributed_fft3d.cc)
target_link_libraries(distributed_fft3d ${PARALLEL_FFTW_LIBS})

add_subdirectory(tests)
