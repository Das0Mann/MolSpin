# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /users/student/pegu0756/molspin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /users/student/pegu0756/molspin

# Include any dependencies generated for this target.
include CMakeFiles/molspin.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/molspin.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/molspin.dir/flags.make

CMakeFiles/molspin.dir/main.cpp.o: CMakeFiles/molspin.dir/flags.make
CMakeFiles/molspin.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/student/pegu0756/molspin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/molspin.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/molspin.dir/main.cpp.o -c /users/student/pegu0756/molspin/main.cpp

CMakeFiles/molspin.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/molspin.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/student/pegu0756/molspin/main.cpp > CMakeFiles/molspin.dir/main.cpp.i

CMakeFiles/molspin.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/molspin.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/student/pegu0756/molspin/main.cpp -o CMakeFiles/molspin.dir/main.cpp.s

# Object files for target molspin
molspin_OBJECTS = \
"CMakeFiles/molspin.dir/main.cpp.o"

# External object files for target molspin
molspin_EXTERNAL_OBJECTS =

molspin: CMakeFiles/molspin.dir/main.cpp.o
molspin: CMakeFiles/molspin.dir/build.make
molspin: libmolspin_core.a
molspin: CMakeFiles/molspin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/users/student/pegu0756/molspin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable molspin"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/molspin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/molspin.dir/build: molspin

.PHONY : CMakeFiles/molspin.dir/build

CMakeFiles/molspin.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/molspin.dir/cmake_clean.cmake
.PHONY : CMakeFiles/molspin.dir/clean

CMakeFiles/molspin.dir/depend:
	cd /users/student/pegu0756/molspin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/student/pegu0756/molspin /users/student/pegu0756/molspin /users/student/pegu0756/molspin /users/student/pegu0756/molspin /users/student/pegu0756/molspin/CMakeFiles/molspin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/molspin.dir/depend
