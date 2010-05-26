# - FindMPI4PY
# Find mpi4py includes and library
# This module defines:
# MPI4PY_INCLUDE_DIR, where to find mpi4py.h, etc.
# MPI4PY_FOUND

if(NOT MPI4PY_INCLUDE_DIR)
    execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "import mpi4py; print mpi4py.get_include(),"
      OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR_OUTPUT
      RESULT_VARIABLE MPI4PY_COMMAND_RESULT)
    if(MPI4PY_COMMAND_RESULT)
        message("jfa: mpi4py not found")
        set(MPI4PY_FOUND FALSE)
    else(MPI4PY_COMMAND_RESULT)
        string(REPLACE "\n" "" MPI4PY_INCLUDE_DIR ${MPI4PY_INCLUDE_DIR_OUTPUT})
        if (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
            message("jfa: mpi4py matches traceback")
            ## Did not successfully include MPI4PY
            set(MPI4PY_FOUND FALSE)
        else (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
            ## successful
            set(MPI4PY_FOUND TRUE)
            set(MPI4PY_INCLUDE_DIR ${MPI4PY_INCLUDE_DIR} CACHE STRING "mpi4py include path")
        endif (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
    endif(MPI4PY_COMMAND_RESULT)
else(NOT MPI4PY_INCLUDE_DIR)
    set(MPI4PY_FOUND TRUE)
endif(NOT MPI4PY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI4PY DEFAULT_MSG MPI4PY_INCLUDE_DIR)

