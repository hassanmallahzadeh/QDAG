# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.18.4/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.18.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator

# Include any dependencies generated for this target.
include CMakeFiles/quantum_simulator.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/quantum_simulator.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/quantum_simulator.dir/flags.make

CMakeFiles/quantum_simulator.dir/src/main.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/quantum_simulator.dir/src/main.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/main.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/main.cpp

CMakeFiles/quantum_simulator.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/main.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/main.cpp > CMakeFiles/quantum_simulator.dir/src/main.cpp.i

CMakeFiles/quantum_simulator.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/main.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/main.cpp -o CMakeFiles/quantum_simulator.dir/src/main.cpp.s

CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o: src/QMDDpackage.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDpackage.cpp

CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDpackage.cpp > CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.i

CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDpackage.cpp -o CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.s

CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o: src/textFileUtilities.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/textFileUtilities.cpp

CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/textFileUtilities.cpp > CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.i

CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/textFileUtilities.cpp -o CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.s

CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o: src/qcost.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/qcost.cpp

CMakeFiles/quantum_simulator.dir/src/qcost.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/qcost.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/qcost.cpp > CMakeFiles/quantum_simulator.dir/src/qcost.cpp.i

CMakeFiles/quantum_simulator.dir/src/qcost.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/qcost.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/qcost.cpp -o CMakeFiles/quantum_simulator.dir/src/qcost.cpp.s

CMakeFiles/quantum_simulator.dir/src/timing.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/timing.cpp.o: src/timing.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/quantum_simulator.dir/src/timing.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/timing.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/timing.cpp

CMakeFiles/quantum_simulator.dir/src/timing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/timing.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/timing.cpp > CMakeFiles/quantum_simulator.dir/src/timing.cpp.i

CMakeFiles/quantum_simulator.dir/src/timing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/timing.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/timing.cpp -o CMakeFiles/quantum_simulator.dir/src/timing.cpp.s

CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o: src/QMDDcircuit.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcircuit.cpp

CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcircuit.cpp > CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.i

CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcircuit.cpp -o CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.s

CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o: src/QMDDcomplexD.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcomplexD.cpp

CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcomplexD.cpp > CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.i

CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDcomplexD.cpp -o CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.s

CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o: CMakeFiles/quantum_simulator.dir/flags.make
CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o: src/QMDDreorder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o -c /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDreorder.cpp

CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDreorder.cpp > CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.i

CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/src/QMDDreorder.cpp -o CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.s

# Object files for target quantum_simulator
quantum_simulator_OBJECTS = \
"CMakeFiles/quantum_simulator.dir/src/main.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/timing.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o" \
"CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o"

# External object files for target quantum_simulator
quantum_simulator_EXTERNAL_OBJECTS =

quantum_simulator: CMakeFiles/quantum_simulator.dir/src/main.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/QMDDpackage.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/textFileUtilities.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/qcost.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/timing.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/QMDDcircuit.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/QMDDcomplexD.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/src/QMDDreorder.cpp.o
quantum_simulator: CMakeFiles/quantum_simulator.dir/build.make
quantum_simulator: CMakeFiles/quantum_simulator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable quantum_simulator"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/quantum_simulator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/quantum_simulator.dir/build: quantum_simulator

.PHONY : CMakeFiles/quantum_simulator.dir/build

CMakeFiles/quantum_simulator.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/quantum_simulator.dir/cmake_clean.cmake
.PHONY : CMakeFiles/quantum_simulator.dir/clean

CMakeFiles/quantum_simulator.dir/depend:
	cd /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator /Users/hassanmallahzadeh/UBCLife/CPSC448/QDAG/AncientAustrianSimulator/CMakeFiles/quantum_simulator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/quantum_simulator.dir/depend

