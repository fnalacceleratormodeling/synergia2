macro(add_mpi_test tname np)

    if (${Kokkos_ENABLE_OPENMP})
        add_test(NAME ${tname}_${np}
            COMMAND ${CMAKE_COMMAND} -E env 
                OMP_NUM_THREADS=1
                OMP_PROC_BIND=spread
                OMP_PLACES=threads
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${tname}> ${MPIEXEC_POSTFLAGS})
    else()
        add_test(NAME ${tname}_${np}
            COMMAND ${CMAKE_COMMAND} -E env
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${tname}> ${MPIEXEC_POSTFLAGS})
    endif()

endmacro()
