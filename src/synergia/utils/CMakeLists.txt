add_library(synergia_parallel_utils 
    parallel_utils.cc 
    commxx.cc 
    #commxx_divider.cc 
    logger.cc
    simple_timer.cc
    base64.cpp
    )
target_link_libraries(synergia_parallel_utils 
    #${Boost_SERIALIZATION_LIBRARIES}
    )

add_library(synergia_test_main STATIC catch_test_main.cc)
target_link_libraries(synergia_test_main kokkos)

add_library(synergia_serialization 
    #serialization_files.cc
    cereal_files.cc
    )

target_link_libraries(synergia_serialization 
#    ${Boost_FILESYSTEM_LIBRARIES}
#    ${Boost_SYSTEM_LIBRARIES} 
#    ${Boost_SERIALIZATION_LIBRARY} 
    synergia_parallel_utils
    )

add_library(synergia_hdf5_utils 
    #hdf5_writer.cc 
    hdf5_serial_writer.cc
    #hdf5_chunked_array2d_writer.cc 
    hdf5_file.cc
    )
target_link_libraries(synergia_hdf5_utils 
    ${HDF5_LIBRARIES} 
    #synergia_serialization
    synergia_parallel_utils
    kokkos
    )

if (BUILD_PYTHON_BINDINGS)
    #add_python_extension(parallel_utils parallel_utils_wrap.cc)
    #target_link_libraries(parallel_utils synergia_parallel_utils synergia_hdf5_utils ${Boost_LIBRARIES})

    #add_python_extension(pylsexpr lsexpr_wrap.cc)
    #target_link_libraries(pylsexpr lsexpr ${Boost_LIBRARIES})

    pybind11_add_module(parallel_utils parallel_utils_pywrap.cc)
    target_link_libraries(parallel_utils PRIVATE synergia_parallel_utils)

    pybind11_add_module(pylsexpr lsexpr_pywrap.cc)
    target_link_libraries(pylsexpr PRIVATE lsexpr)

    copy_file(__init__.py parallel_utils)
endif ()

if(${KOKKOS_ENABLE_CUDA})
    set(FFT_SRC distributed_fft2d_cuda.cc)
    set(FFT_LIB cufft)
else()
    set(FFT_SRC distributed_fft2d_fftw.cc)
    set(FFT_LIB ${PARALLEL_FFTW_LIBRARIES})
endif()

add_library(synergia_distributed_fft 
    ${FFT_SRC}
    )

target_link_libraries(synergia_distributed_fft 
    ${FFT_LIB}
    #synergia_serialization 
    synergia_parallel_utils
    )

add_library(synergia_command_line command_line_arg.cc)

add_library(lsexpr lsexpr.cc)

install(TARGETS
    synergia_parallel_utils
    synergia_hdf5_utils
    synergia_distributed_fft
    #synergia_serialization
    synergia_command_line
    lsexpr
    DESTINATION ${LIB_INSTALL_DIR})

install(FILES
    boost_test_mpi_fixture.h
    comm_converter.h
    command_line_arg.h
    commxx.h
    commxx_divider.h
    container_conversions.h
    distributed_fft3d.h
    distributed_fft2d.h
    fast_int_floor.h
    floating_point.h
    gsvector.h
    hdf5_chunked_array2d_writer.h
    hdf5_misc.h
    hdf5_file.h
    hdf5_file.tcc
    hdf5_serial_writer.h
    hdf5_serial_writer.tcc
    hdf5_writer.h
    hdf5_writer.tcc
    complex_error_function.h
    multi_array_assert.h
    multi_array_check_equal.h
    multi_array_offsets.h
    multi_array_print.h
    multi_array_serialization.h
    multi_array_to_string.h
    multi_array_typedefs.h
    numpy_multi_array_converter.h
    numpy_multi_ref_converter.h
    parallel_utils.h
    simple_timer.h
    serialization.h
    serialization_files.h
    cereal.h
    cereal_files.h
    digits.h
    logger.h
    lsexpr.h
    synergia_omp.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/utils)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
    install(TARGETS
        parallel_utils
        #pylsexpr
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/utils)
endif ()

#add_subdirectory(tests)

