macro(add_mpi_test tname np)

    add_test(NAME ${tname}_${np}
        COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} $<TARGET_FILE:${tname}>)

endmacro()
