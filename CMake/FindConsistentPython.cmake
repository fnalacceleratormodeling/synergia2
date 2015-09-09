# FindConsistentPython
#
# This module defines:
#    PYTHON_EXECUTABLE
#    PYTHON_INCLUDE_DIR
#    PYTHON_LIB_DIR
#    PYTHON_LIBRARY
# Each variable may be overridden by the end user.
# The module will attempt to find values for
# unset variables that are consistent with the set ones.

find_package(PythonInterp REQUIRED)
find_package(PythonInclude REQUIRED)

if(NOT PYTHON_LIBRARY)
    find_package(PythonLibDir REQUIRED)
    set(PYTHON_CURRENT_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
    set(PYTHON_CURRENT_VERSION_NO_DOTS "${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
    find_library(PYTHON_LIBRARY
        NAMES python${PYTHON_CURRENT_VERSION_NO_DOTS} python${PYTHON_CURRENT_VERSION}
        NO_SYSTEM_ENVIRONMENT_PATH
        HINTS ${PYTHON_LIB_DIR}
    )
    message(STATUS "FoundPythonLibrary: ${PYTHON_LIBRARY}")
endif(NOT PYTHON_LIBRARY)
