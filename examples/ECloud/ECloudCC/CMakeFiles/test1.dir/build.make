# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.4

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC

# Include any dependencies generated for this target.
include CMakeFiles/test1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/test1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test1.dir/flags.make

CMakeFiles/test1.dir/depend.make.mark: CMakeFiles/test1.dir/flags.make
CMakeFiles/test1.dir/depend.make.mark: test1.cc

CMakeFiles/test1.dir/test1.o: CMakeFiles/test1.dir/flags.make
CMakeFiles/test1.dir/test1.o: test1.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/test1.dir/test1.o"
	/usr/bin/c++   $(CXX_FLAGS) -o CMakeFiles/test1.dir/test1.o -c /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC/test1.cc

CMakeFiles/test1.dir/test1.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test1.dir/test1.i"
	/usr/bin/c++  $(CXX_FLAGS) -E /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC/test1.cc > CMakeFiles/test1.dir/test1.i

CMakeFiles/test1.dir/test1.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test1.dir/test1.s"
	/usr/bin/c++  $(CXX_FLAGS) -S /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC/test1.cc -o CMakeFiles/test1.dir/test1.s

CMakeFiles/test1.dir/test1.o.requires:

CMakeFiles/test1.dir/test1.o.provides: CMakeFiles/test1.dir/test1.o.requires
	$(MAKE) -f CMakeFiles/test1.dir/build.make CMakeFiles/test1.dir/test1.o.provides.build

CMakeFiles/test1.dir/test1.o.provides.build: CMakeFiles/test1.dir/test1.o

CMakeFiles/test1.dir/depend: CMakeFiles/test1.dir/depend.make.mark

CMakeFiles/test1.dir/depend.make.mark:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --magenta --bold "Scanning dependencies of target test1"
	cd /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC /local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC/CMakeFiles/test1.dir/DependInfo.cmake

# Object files for target test1
test1_OBJECTS = \
"CMakeFiles/test1.dir/test1.o"

# External object files for target test1
test1_EXTERNAL_OBJECTS =

test1: CMakeFiles/test1.dir/test1.o
test1: libECloudLib.a
test1: CMakeFiles/test1.dir/build.make
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable test1"
	$(CMAKE_COMMAND) -P CMakeFiles/test1.dir/cmake_clean_target.cmake
	/usr/bin/c++      -fPIC $(test1_OBJECTS) $(test1_EXTERNAL_OBJECTS)  -o test1 -rdynamic -L/local/lebrun/Synergia/cca/install/lib/libpython2.4 -L/local/lebrun/Synergia/cca/install/lib -L/local/lebrun/Tech-Xlib/install/lib -L/local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC -Wl,-Bstatic -lECloudLib -Wl,-Bdynamic -L/local/lebrun/Synergia/cca/install/lib -lgsl -lgslcblas -lm -ltxionpack -ltxegenelec -ltxrand -lboost_python-gcc -lpython2.4 -ltk -ltcl -Wl,-rpath,/local/lebrun/Synergia/cca/install/lib/libpython2.4:/local/lebrun/Synergia/cca/install/lib:/local/lebrun/Tech-Xlib/install/lib:/local/lebrun/Synergia/cca/build/synergia2/examples/ECloud/ECloudCC 

# Rule to build all files generated by this target.
CMakeFiles/test1.dir/build: test1

CMakeFiles/test1.dir/requires: CMakeFiles/test1.dir/test1.o.requires

CMakeFiles/test1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test1.dir/cmake_clean.cmake

