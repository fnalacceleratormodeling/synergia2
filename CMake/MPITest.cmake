# OMPI flag for oversubcription: https://github.com/open-mpi/ompi/issues/8955

macro(add_mpi_test tname np)

    if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "OpenMP")
        add_test(NAME ${tname}_${np}
            COMMAND ${CMAKE_COMMAND} -E env 
                OMP_NUM_THREADS=1
                OMP_PROC_BIND=spread
                OMP_PLACES=cores
                OMPI_MCA_rmaps_base_oversubscribe=true
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${tname}> ${MPIEXEC_POSTFLAGS})
    else()
        add_test(NAME ${tname}_${np}
            COMMAND ${CMAKE_COMMAND} -E env
                OMPI_MCA_rmaps_base_oversubscribe=true
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${tname}> ${MPIEXEC_POSTFLAGS})
    endif()

endmacro()
