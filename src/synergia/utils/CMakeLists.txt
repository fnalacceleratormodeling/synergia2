add_library(synergia_parallel_utils parallel_utils.cc commxx.cc)
if (BUILD_PYTHON_BINDINGS)
    add_python_extension(parallel_utils parallel_utils_wrap.cc)
   target_link_libraries(parallel_utils synergia_parallel_utils
        ${Boost_LIBRARIES})
endif ()

add_library(synergia_hdf5_utils hdf5_writer.cc hdf5_serial_writer.cc
    hdf5_chunked_array2d_writer.cc hdf5_file.cc)
target_link_libraries(synergia_hdf5_utils ${HDF5_LIBRARIES})

add_library(synergia_distributed_fft distributed_fft3d.cc
    distributed_fft2d.cc)
target_link_libraries(synergia_distributed_fft ${PARALLEL_FFTW_LIBRARIES} synergia_parallel_utils)

add_library(synergia_command_line command_line_arg.cc)

install(TARGETS
    synergia_parallel_utils
    synergia_hdf5_utils
    synergia_distributed_fft
    DESTINATION lib)
install(FILES
    boost_test_mpi_fixture.h
    comm_converter.h
    command_line_arg.h
    commxx.h
    container_conversions.h
    distributed_fft3d.h
    distributed_fft2d.h
    fast_int_floor.h
    floating_point.h
    hdf5_chunked_array2d_writer.h
    hdf5_misc.h
    hdf5_serial_writer.h
    hdf5_serial_writer.tcc
    hdf5_writer.h
    hdf5_writer.tcc
    complex_error_function.h
    multi_array_check_equal.h
    multi_array_offsets.h
    multi_array_print.h
    multi_array_serialization.h
    multi_array_typedefs.h
    numpy_multi_ref_converter.h
    parallel_utils.h
    simple_timer.h
    xml_serialization.h
    DESTINATION include/synergia/utils)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/utils)
    install(TARGETS
        parallel_utils
        DESTINATION lib/synergia/utils)
endif ()

add_subdirectory(tests)
