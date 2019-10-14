# - FindMPI4PY
# Find mpi4py includes and library
# This module defines:
# MPI4PY_INCLUDE_DIR, where to find mpi4py.h, etc.
# MPI4PY_FOUND

if(NOT MPI4PY_INCLUDE_DIR)
#     if(EXTRA_PYTHONPATH)
#         set(OLD_PYTHONPATH $ENV{PYTHONPATH})
#         set(ENV{PYTHONPATH} "${EXTRA_PYTHONPATH}:$ENV{PYTHONPATH}")
#     endif(EXTRA_PYTHONPATH)
    execute_process(COMMAND "${MY_PYTHON_EXECUTABLE}"
      "-c" "from __future__ import print_function; import mpi4py; print(mpi4py.get_include())"
      OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
      RESULT_VARIABLE MPI4PY_NOT_FOUND
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(MPI4PY_NOT_FOUND)
      set(INTERNAL_MPI4PY_FOUND FALSE)
        message("FindMPI4PY: mpi4py not found")
        set(MPI4PY_FOUND FALSE)
    else()
        if (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
            message("FindMPI4PY: mpi4py matches traceback")
            ## Did not successfully include MPI4PY
            set(MPI4PY_FOUND FALSE)
        else (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
            ## successful
            set(INTERNALMPI4PY_FOUND TRUE)
            set(MPI4PY_INCLUDE_DIR ${MPI4PY_INCLUDE_DIR} CACHE STRING "mpi4py include path")
        endif (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
    endif()
    set(MPI4PY_FOUND ${INTERNAL_MPI4PY_FOUND} CACHE BOOL "mpi4py found" FORCE)
#     if(EXTRA_PYTHONPATH)
#         set(ENV{PYTHONPATH} ${OLD_PYTHONPATH})
#     endif(EXTRA_PYTHONPATH)
else(NOT MPI4PY_INCLUDE_DIR)
    set(MPI4PY_FOUND TRUE)
endif(NOT MPI4PY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI4PY DEFAULT_MSG MPI4PY_INCLUDE_DIR)

