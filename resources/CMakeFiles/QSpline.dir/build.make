# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/czhou/Projects/QuinticSpline

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/czhou/Projects/QuinticSpline/resources

# Include any dependencies generated for this target.
include CMakeFiles/QSpline.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/QSpline.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/QSpline.dir/flags.make

CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o: CMakeFiles/QSpline.dir/flags.make
CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o: ../src/quinticSpline.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/czhou/Projects/QuinticSpline/resources/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o -c /home/czhou/Projects/QuinticSpline/src/quinticSpline.cpp

CMakeFiles/QSpline.dir/src/quinticSpline.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/QSpline.dir/src/quinticSpline.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/czhou/Projects/QuinticSpline/src/quinticSpline.cpp > CMakeFiles/QSpline.dir/src/quinticSpline.cpp.i

CMakeFiles/QSpline.dir/src/quinticSpline.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/QSpline.dir/src/quinticSpline.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/czhou/Projects/QuinticSpline/src/quinticSpline.cpp -o CMakeFiles/QSpline.dir/src/quinticSpline.cpp.s

# Object files for target QSpline
QSpline_OBJECTS = \
"CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o"

# External object files for target QSpline
QSpline_EXTERNAL_OBJECTS =

libQSpline.a: CMakeFiles/QSpline.dir/src/quinticSpline.cpp.o
libQSpline.a: CMakeFiles/QSpline.dir/build.make
libQSpline.a: CMakeFiles/QSpline.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/czhou/Projects/QuinticSpline/resources/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libQSpline.a"
	$(CMAKE_COMMAND) -P CMakeFiles/QSpline.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/QSpline.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/QSpline.dir/build: libQSpline.a

.PHONY : CMakeFiles/QSpline.dir/build

CMakeFiles/QSpline.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/QSpline.dir/cmake_clean.cmake
.PHONY : CMakeFiles/QSpline.dir/clean

CMakeFiles/QSpline.dir/depend:
	cd /home/czhou/Projects/QuinticSpline/resources && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/czhou/Projects/QuinticSpline /home/czhou/Projects/QuinticSpline /home/czhou/Projects/QuinticSpline/resources /home/czhou/Projects/QuinticSpline/resources /home/czhou/Projects/QuinticSpline/resources/CMakeFiles/QSpline.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/QSpline.dir/depend

