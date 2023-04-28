macro(add_py_test tname)

  add_custom_target(pytest_${tname})

  add_test(
    NAME ${tname}
    COMMAND
      ${CMAKE_COMMAND} -E env OMPI_MCA_rmaps_base_oversubscribe=true
      PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${Python_EXECUTABLE} -m
      pytest -vv ${CMAKE_CURRENT_BINARY_DIR}/${tname} ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${SYNERGIA2_BINARY_DIR}/src)

  copy_file(${tname} pytest_${tname})

endmacro()
