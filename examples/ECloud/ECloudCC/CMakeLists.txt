project(ECloudCC)

# This should be set from above.. 
set(SYNERGIA2_SOURCE_DIR /local/lebrun/Synergia/cca/build/synergia2/)

enable_testing()

INCLUDE(FindPythonInterp)
INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})
message(" Python libraries are " ${PYTHON_LIBRARIES})
INCLUDE(AddPythonExtension.cmake)

FIND_PACKAGE(Boost)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
message(" Boost libraries are " ${Boost_LIBRARIES})
set(BOOST_LIBS -lboost_python-gcc)
#
# Chef stuff is needed... should work, but it does not... 
#
#INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindCHEF.cmake)
set (CHEF_INCLUDE_DIRS /local/lebrun/Synergia/cca/install/include/)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
set (CHEF_LIBRARY_DIRS /local/lebrun/Synergia/cca/install/lib/)
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})
#

include(FindGSL.cmake)

# Add the Tech-X stuff,explicit directories. 
# ( to be included in contract later on) 

set (TECHX_INCLUDE_DIRS /local/lebrun/Tech-Xlib/install/include)
INCLUDE_DIRECTORIES(${TECHX_INCLUDE_DIRS})
set (TECHX_LIBRARY_DIRS /local/lebrun/Tech-Xlib/install/lib)
LINK_DIRECTORIES(${TECHX_LIBRARY_DIRS})
#
# Add the Synergia stuff, if invoked from above, this should not be needed. 
#
INCLUDE_DIRECTORIES(${SYNERGIA2_SOURCE_DIR}/s2_fish)
#
# Python Module
#
add_python_extension(ECloudPy ECloudPy.cc util1.cc RKIntegrators.cc 
                     RKIntegratorWithElec.cc MIDipoleEdge.cc VorpalYeeField.cc)
target_link_libraries(ECloudPy -ltxionpack -ltxegenelec -ltxrand -lhdf5_cpp -lhdf5)
target_link_libraries(ECloudPy ${BOOST_LIBS} -lpython2.4)
target_link_libraries(ECloudPy ${CHEF_LIBRARY_DIRS} libbasic_toolkit.so libfftw.so)
target_link_libraries(ECloudPy ${GSL_CONFIG_LIBS})
target_link_libraries(ECloudPy -ltk -ltcl)

add_library(ECloudLib util1.cc)
target_link_libraries(ECloudLib ${GSL_CONFIG_LIBS} ${CHEF_LIBRARY_DIRS})
#
add_executable(test1 test1.cc)
target_link_libraries(test1 ECloudLib)
target_link_libraries(ECloudLib -ltxionpack -ltxegenelec -ltxrand)
target_link_libraries(ECloudLib ${BOOST_LIBS} -lpython2.4.so)
target_link_libraries(ECloudLib -lbasic_toolkit -lfftw -lrfftw -lfftw3)
target_link_libraries(ECloudLib -ltk -ltcl)
#
add_executable(testFieldIntGSL testFieldIntGSL.cc RKIntegrators.cc MIDipoleEdge.cc util1.cc VorpalYeeField.cc)
target_link_libraries(testFieldIntGSL ${GSL_CONFIG_LIBS})
target_link_libraries(testFieldIntGSL -lbasic_toolkit -lfftw -lrfftw -lfftw3)
target_link_libraries(testFieldIntGSL -ltxionpack -ltxegenelec -ltxrand -lhdf5_cpp -lhdf5)
target_link_libraries(testFieldIntGSL ${BOOST_LIBS} -lpython2.4)
target_link_libraries(testFieldIntGSL ${GSL_CONFIG_LIBS})
target_link_libraries(testFieldIntGSL -ltk -ltcl)
#
add_executable(test2 test2.cc biGaussFunc.cc)
target_link_libraries(test2 ${GSL_CONFIG_LIBS})
#
add_executable(test3 test3.cc VorpalYeeField.cc)
target_link_libraries(test3 -lhdf5_cpp -lhdf5 )
#
add_executable(test4 test4.cc RKIntegrators.cc MIDipoleEdge.cc util1.cc VorpalYeeField.cc VorpalElectrons.cc)
target_link_libraries(test4 ${GSL_CONFIG_LIBS})
target_link_libraries(test4 -lbasic_toolkit -lfftw -lrfftw -lfftw3)
target_link_libraries(test4 -ltxionpack -ltxegenelec -ltxrand -lhdf5_cpp -lhdf5)
target_link_libraries(test4 ${BOOST_LIBS} -lpython2.4)
target_link_libraries(test4 ${GSL_CONFIG_LIBS})
target_link_libraries(test4 -lhdf5_cpp -lhdf5 )
target_link_libraries(test4 -ltk -ltcl)
#
add_executable(testBeamKick testBeamKick.cc RKIntegrators.cc MIDipoleEdge.cc util1.cc VorpalYeeField.cc)
target_link_libraries(testBeamKick ${})
target_link_libraries(testBeamKick -lbasic_toolkit -lfftw -lrfftw -lfftw3)
target_link_libraries(testBeamKick -ltxionpack -ltxegenelec -ltxrand)
target_link_libraries(testBeamKick ${GSL_CONFIG_LIBS})
target_link_libraries(testBeamKick ${BOOST_LIBS} -lpython2.4)
target_link_libraries(testBeamKick ${GSL_CONFIG_LIBS})
target_link_libraries(testBeamKick -ltk -ltcl)
#
