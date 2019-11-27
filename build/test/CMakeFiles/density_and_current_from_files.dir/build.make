# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/jgarcia/codes/wannier2sparse

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/jgarcia/codes/wannier2sparse/build

# Include any dependencies generated for this target.
include test/CMakeFiles/density_and_current_from_files.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/density_and_current_from_files.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/density_and_current_from_files.dir/flags.make

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o: test/CMakeFiles/density_and_current_from_files.dir/flags.make
test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o: ../test/density_and_current_from_files.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/jgarcia/codes/wannier2sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o"
	cd /data/jgarcia/codes/wannier2sparse/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o -c /data/jgarcia/codes/wannier2sparse/test/density_and_current_from_files.cpp

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.i"
	cd /data/jgarcia/codes/wannier2sparse/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/jgarcia/codes/wannier2sparse/test/density_and_current_from_files.cpp > CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.i

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.s"
	cd /data/jgarcia/codes/wannier2sparse/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/jgarcia/codes/wannier2sparse/test/density_and_current_from_files.cpp -o CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.s

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.requires:

.PHONY : test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.requires

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.provides: test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.requires
	$(MAKE) -f test/CMakeFiles/density_and_current_from_files.dir/build.make test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.provides.build
.PHONY : test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.provides

test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.provides.build: test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o


# Object files for target density_and_current_from_files
density_and_current_from_files_OBJECTS = \
"CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o"

# External object files for target density_and_current_from_files
density_and_current_from_files_EXTERNAL_OBJECTS =

test/density_and_current_from_files/density_and_current_from_files: test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o
test/density_and_current_from_files/density_and_current_from_files: test/CMakeFiles/density_and_current_from_files.dir/build.make
test/density_and_current_from_files/density_and_current_from_files: lib/libwannierlib.a
test/density_and_current_from_files/density_and_current_from_files: test/CMakeFiles/density_and_current_from_files.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/jgarcia/codes/wannier2sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable density_and_current_from_files/density_and_current_from_files"
	cd /data/jgarcia/codes/wannier2sparse/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/density_and_current_from_files.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/density_and_current_from_files.dir/build: test/density_and_current_from_files/density_and_current_from_files

.PHONY : test/CMakeFiles/density_and_current_from_files.dir/build

test/CMakeFiles/density_and_current_from_files.dir/requires: test/CMakeFiles/density_and_current_from_files.dir/density_and_current_from_files.cpp.o.requires

.PHONY : test/CMakeFiles/density_and_current_from_files.dir/requires

test/CMakeFiles/density_and_current_from_files.dir/clean:
	cd /data/jgarcia/codes/wannier2sparse/build/test && $(CMAKE_COMMAND) -P CMakeFiles/density_and_current_from_files.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/density_and_current_from_files.dir/clean

test/CMakeFiles/density_and_current_from_files.dir/depend:
	cd /data/jgarcia/codes/wannier2sparse/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/jgarcia/codes/wannier2sparse /data/jgarcia/codes/wannier2sparse/test /data/jgarcia/codes/wannier2sparse/build /data/jgarcia/codes/wannier2sparse/build/test /data/jgarcia/codes/wannier2sparse/build/test/CMakeFiles/density_and_current_from_files.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/density_and_current_from_files.dir/depend

