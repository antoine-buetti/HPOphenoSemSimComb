# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build

# Include any dependencies generated for this target.
include test/CMakeFiles/HPOphenoSemSimCombTest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/HPOphenoSemSimCombTest.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make

test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/findMICATest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/findMICATest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/findMICATest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/findMICATest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/frameDataTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/frameDataTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/frameDataTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/frameDataTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/grepHPOforCorrespondingGenesTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/grepHPOforCorrespondingGenesTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/grepHPOforCorrespondingGenesTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/grepHPOforCorrespondingGenesTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/main.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/main.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/main.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/main.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetAssociatedHPOtermsTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetAssociatedHPOtermsTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetAssociatedHPOtermsTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetAssociatedHPOtermsTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetChildrenTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetChildrenTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetChildrenTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetChildrenTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetParentsTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetParentsTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetParentsTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/queryNodeGetParentsTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/semSimAndICTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/semSimAndICTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/semSimAndICTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/semSimAndICTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkDownTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkDownTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkDownTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkDownTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkUpTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkUpTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkUpTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/treeWalkUpTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.s

test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/flags.make
test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o: /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/uniqueNodeAndAllDescendantsTest.cpp
test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o: test/CMakeFiles/HPOphenoSemSimCombTest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o -MF CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o.d -o CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o -c /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/uniqueNodeAndAllDescendantsTest.cpp

test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.i"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/uniqueNodeAndAllDescendantsTest.cpp > CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.i

test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.s"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test/uniqueNodeAndAllDescendantsTest.cpp -o CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.s

# Object files for target HPOphenoSemSimCombTest
HPOphenoSemSimCombTest_OBJECTS = \
"CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o" \
"CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o"

# External object files for target HPOphenoSemSimCombTest
HPOphenoSemSimCombTest_EXTERNAL_OBJECTS =

test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/findMICATest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/frameDataTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/grepHPOforCorrespondingGenesTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/main.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetAssociatedHPOtermsTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryGeneNameGetOutputSemSimDiseaseIDTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetChildrenTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/queryNodeGetParentsTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/semSimAndICTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkDownTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/treeWalkUpTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/uniqueNodeAndAllDescendantsTest.cpp.o
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/build.make
test/HPOphenoSemSimCombTest: src/libHPOphenoSemSimComb_lib.a
test/HPOphenoSemSimCombTest: lib/libgtest.a
test/HPOphenoSemSimCombTest: test/CMakeFiles/HPOphenoSemSimCombTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable HPOphenoSemSimCombTest"
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HPOphenoSemSimCombTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/HPOphenoSemSimCombTest.dir/build: test/HPOphenoSemSimCombTest
.PHONY : test/CMakeFiles/HPOphenoSemSimCombTest.dir/build

test/CMakeFiles/HPOphenoSemSimCombTest.dir/clean:
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test && $(CMAKE_COMMAND) -P CMakeFiles/HPOphenoSemSimCombTest.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/HPOphenoSemSimCombTest.dir/clean

test/CMakeFiles/HPOphenoSemSimCombTest.dir/depend:
	cd /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/test /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test /home/runner/work/HPOphenoSemSimComb/HPOphenoSemSimComb/build/test/CMakeFiles/HPOphenoSemSimCombTest.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/CMakeFiles/HPOphenoSemSimCombTest.dir/depend

