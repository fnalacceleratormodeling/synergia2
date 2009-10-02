# - Find mpi4py

# Find the mpi4py includes
# This module defines
#  MPI4PY_INCLUDE_DIR, where to find MPI4PY/arrayobject.h, etc.
#  MPI4PY_FOUND, If false, do not try to use MPI4PY headers.

#if (MPI4PY_INCLUDE_DIR)
  # in cache already
#  set (MPI4PY_FIND_QUIETLY TRUE)
#endif (MPI4PY_INCLUDE_DIR)

EXEC_PROGRAM ("${PYTHON_EXECUTABLE}"
  ARGS "-c 'import mpi4py; print mpi4py.get_include()'"
  OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
  RETURN_VALUE MPI4PY_NOT_FOUND)

if (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
# Did not successfully include MPI4PY
  set(MPI4PY_FOUND FALSE)
else (MPI4PY_INCLUDE_DIR MATCHES "Traceback")
# successful
  set (MPI4PY_FOUND TRUE)
  set (MPI4PY_INCLUDE_DIR ${MPI4PY_INCLUDE_DIR} CACHE STRING "mpi4py include path")
endif (MPI4PY_INCLUDE_DIR MATCHES "Traceback")

if (MPI4PY_FOUND)
  if (NOT MPI4PY_FIND_QUIETLY)
    message (STATUS "mpi4py headers found")
  endif (NOT MPI4PY_FIND_QUIETLY)
else (MPI4PY_FOUND)
  if (MPI4PY_FIND_REQUIRED)
    message (FATAL_ERROR "mpi4py headers missing")
  endif (MPI4PY_FIND_REQUIRED)
endif (MPI4PY_FOUND)

MARK_AS_ADVANCED (MPI4PY_INCLUDE_DIR)
