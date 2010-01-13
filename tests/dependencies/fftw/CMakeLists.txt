add_python_extension(trivial_fftw trivial_fftw.cc)
target_link_libraries(trivial_fftw ${Boost_LIBRARIES} 
    ${FFTW2_R_MPI_LIB} ${FFTW2_MPI_LIB} ${FFTW2_R_LIB} ${FFTW2_LIB})

