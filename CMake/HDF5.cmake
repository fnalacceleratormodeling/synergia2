message(STATUS "fetching hdf5")
set(HDF5_ENABLE_PARALLEL ON)
set(BUILD_SHARED_LIBS ON)
set(HDF5_EXTERNALLY_CONFIGURED 1)
include(FetchContent)
FetchContent_Declare(
  hdf5
  URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/CMake-hdf5-1.12.2.tar.gz)
  #GIT_REPOSITORY https://github.com/HDFGroup/hdf5
  #GIT_TAG "hdf5-1_12_2")
FetchContent_MakeAvailable(hdf5)
#add_library(hdf5 INTERFACE ${hdf5_SOURCE_DIR})
#add_library(hdf5::hdf5 INTERFACE IMPORTED GLOBAL)
#target_include_directories(hdf5::hdf5 INTERFACE "${HDF5_INCLUDE_DIRS}")
#target_link_libraries(hdf5::hdf5 INTERFACE "${HDF5_LIBRARIES}")

#add_dependencies(HDF5::HDF5 HDF5)
#find_package(hdf5 REQUIRED
#	${hdf5_BUILD_DIR})

