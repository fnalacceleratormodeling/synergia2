macro(add_py_test tname)

  add_custom_target(pytest_${tname})

  add_test(
    NAME ${tname}
    COMMAND
      ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
      ${Python_EXECUTABLE} -m pytest -vv ${CMAKE_CURRENT_BINARY_DIR}/${tname}
      ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${SYNERGIA2_BINARY_DIR}/src)

  copy_file(${tname} pytest_${tname})

endmacro()
