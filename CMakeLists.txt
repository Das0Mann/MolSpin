# Set cmake version requirement
cmake_minimum_required(VERSION 3.14)
project(MolSpin)

option(LINK_OPENBLAS "Whether OpenBLAS should be linked." ON)
option(LINK_LAPACK "Whether LAPACK should be linked." ON)
option(ENABLE_OPENMP "Whether OpenMP should be available from the MolSpin code." ON)

# Compiler options
set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wall -O2)
# Paths
set(PATH_SOURCE_MOLSPIN_MAIN "${CMAKE_CURRENT_SOURCE_DIR}")
set(PATH_SOURCE_SPINAPI "${CMAKE_CURRENT_SOURCE_DIR}/SpinAPI")
set(PATH_SOURCE_MSDPARSER "${CMAKE_CURRENT_SOURCE_DIR}/MSDParser")
set(PATH_SOURCE_RUNSECTION "${CMAKE_CURRENT_SOURCE_DIR}/RunSection")
set(PATH_SOURCE_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tests")
# Include modules
include_directories(${PATH_SOURCE_SPINAPI})
include_directories(${PATH_SOURCE_MSDPARSER})
include_directories(${PATH_SOURCE_RUNSECTION})
include_directories("${PATH_SOURCE_RUNSECTION}/Tasks")
include_directories("${PATH_SOURCE_RUNSECTION}/Tasks/Custom")
include_directories("${PATH_SOURCE_RUNSECTION}/Actions")
include_directories(${PATH_SOURCE_TESTS})
# -------------------------------------------------
# SpinAPI source code
# -------------------------------------------------
set(SOURCES_MODULE_SPINAPI
	${PATH_SOURCE_SPINAPI}/Spin.h
	${PATH_SOURCE_SPINAPI}/Spin.cpp
	${PATH_SOURCE_SPINAPI}/SpinSpace.h
	${PATH_SOURCE_SPINAPI}/SpinSpace.cpp
	${PATH_SOURCE_SPINAPI}/SpinSystem.h
	${PATH_SOURCE_SPINAPI}/SpinSystem.cpp
	${PATH_SOURCE_SPINAPI}/Interaction.h
	${PATH_SOURCE_SPINAPI}/Interaction.cpp
	${PATH_SOURCE_SPINAPI}/Transition.h
	${PATH_SOURCE_SPINAPI}/Transition.cpp
	${PATH_SOURCE_SPINAPI}/State.h
	${PATH_SOURCE_SPINAPI}/State.cpp
	${PATH_SOURCE_SPINAPI}/Operator.h
	${PATH_SOURCE_SPINAPI}/Operator.cpp
	${PATH_SOURCE_SPINAPI}/StandardOutput.h
	${PATH_SOURCE_SPINAPI}/StandardOutput.cpp
	${PATH_SOURCE_SPINAPI}/Tensor.h
	${PATH_SOURCE_SPINAPI}/Tensor.cpp
	${PATH_SOURCE_SPINAPI}/Trajectory.h
	${PATH_SOURCE_SPINAPI}/Trajectory.cpp
	${PATH_SOURCE_SPINAPI}/SpinAPIDefines.h
	${PATH_SOURCE_SPINAPI}/SpinAPIfwd.h
)
# -------------------------------------------------
# MSD Parser source code
# -------------------------------------------------
set(SOURCES_MODULE_MSDPARSER
	${PATH_SOURCE_MSDPARSER}/MSDParser.h
	${PATH_SOURCE_MSDPARSER}/MSDParser.cpp
	${PATH_SOURCE_MSDPARSER}/ObjectParser.h
	${PATH_SOURCE_MSDPARSER}/ObjectParser.cpp
	${PATH_SOURCE_MSDPARSER}/FileReader.h
	${PATH_SOURCE_MSDPARSER}/FileReader.cpp
	${PATH_SOURCE_MSDPARSER}/MSDParserDefines.h
	${PATH_SOURCE_MSDPARSER}/MSDParserfwd.h
)
# -------------------------------------------------
# Run Section source code
# -------------------------------------------------
set(SOURCES_MODULE_RUNSECTION
	# Main Run Section functionality
	${PATH_SOURCE_RUNSECTION}/Action.h
	${PATH_SOURCE_RUNSECTION}/Action.cpp
	${PATH_SOURCE_RUNSECTION}/ActionTarget.h
	${PATH_SOURCE_RUNSECTION}/BasicTask.h
	${PATH_SOURCE_RUNSECTION}/BasicTask.cpp
	${PATH_SOURCE_RUNSECTION}/OutputHandler.h
	${PATH_SOURCE_RUNSECTION}/OutputHandler.cpp
	${PATH_SOURCE_RUNSECTION}/RunSection.h
	${PATH_SOURCE_RUNSECTION}/RunSection.cpp
	${PATH_SOURCE_RUNSECTION}/Settings.h
	${PATH_SOURCE_RUNSECTION}/Settings.cpp
	${PATH_SOURCE_RUNSECTION}/RunSectionDefines.h
	${PATH_SOURCE_RUNSECTION}/RunSectionfwd.h
	# Actions
	${PATH_SOURCE_RUNSECTION}/Actions/ActionAddScalar.h
	${PATH_SOURCE_RUNSECTION}/Actions/ActionAddScalar.cpp
	${PATH_SOURCE_RUNSECTION}/Actions/ActionAddVector.h
	${PATH_SOURCE_RUNSECTION}/Actions/ActionAddVector.cpp
	${PATH_SOURCE_RUNSECTION}/Actions/ActionMultiplyScalar.h
	${PATH_SOURCE_RUNSECTION}/Actions/ActionMultiplyScalar.cpp
	${PATH_SOURCE_RUNSECTION}/Actions/ActionRotateVector.h
	${PATH_SOURCE_RUNSECTION}/Actions/ActionRotateVector.cpp
	${PATH_SOURCE_RUNSECTION}/Actions/ActionScaleVector.h
	${PATH_SOURCE_RUNSECTION}/Actions/ActionScaleVector.cpp
)
# Tasks are treated separately so they can be separated out into their own library
set(SOURCES_MODULE_RUNSECTION_TASKS
	# Built-in tasks
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskDynamicHSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskDynamicHSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskGammaCompute.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskGammaCompute.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskHamiltonianEigenvalues.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskHamiltonianEigenvalues.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiDynamicHSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiDynamicHSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiStaticSSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiStaticSSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskPeriodicHSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskPeriodicHSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskPeriodicSSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskPeriodicSSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticHSSymmetricDecay.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticHSSymmetricDecay.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticRPOnlyHSSymDec.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticRPOnlyHSSymDec.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSS.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSS.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSTimeEvo.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSTimeEvo.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSSpectra.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfield.h
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfield.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldSparse.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldSparse.cpp
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldTimeEvo.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldTimeEvo.cpp
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldTimeEvoSparse.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSRedfieldTimeEvoSparse.cpp
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSSpectra.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSSpectra.cpp
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSCIDNP.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskStaticSSCIDNP.cpp
	${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiStaticSSRedfieldTimeEvo.h
        ${PATH_SOURCE_RUNSECTION}/Tasks/TaskMultiStaticSSRedfieldTimeEvo.cpp

	# TODO: Add custom tasks here
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskCheckPMP.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskCheckPMP.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskDensityTimeEvo.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskDensityTimeEvo.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskInformation.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskInformation.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskOCTGPM.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskOCTGPM.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskPMPLambda.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskPMPLambda.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskRPOCT2.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskRPOCT2.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskRPOCT.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskRPOCT.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskStaticSSRedfield.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskStaticSSRedfield.cpp
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskStaticSSRelaxation.h
	#${PATH_SOURCE_RUNSECTION}/Tasks/Custom/TaskStaticSSRelaxation.cpp
)
# -------------------------------------------------
# Sources for targets
# -------------------------------------------------
# Library containing most of the molspin functionality
set(SOURCES_TARGET_LIBRARY_MOLSPIN_CORE ${SOURCES_MODULE_SPINAPI} ${SOURCES_MODULE_MSDPARSER} ${SOURCES_MODULE_RUNSECTION} ${SOURCES_MODULE_RUNSECTION_TASKS})
# The molspin executable
set(SOURCES_TARGET_APPLICATION_MOLSPIN ${PATH_SOURCE_MOLSPIN_MAIN}/main.cpp)
# The test executable
set(SOURCES_TARGET_APPLICATION_TESTS ${PATH_SOURCE_TESTS}/testmain.cpp)
# -------------------------------------------------
# Libraries to link
# -------------------------------------------------
# Core molspin libraries
set(LIBS_MOLSPIN_CORE
	molspin_core
)
# Dependencies
set(LIBS_DEPENDENCIES
	armadillo
)
if(LINK_OPENBLAS)
	list(APPEND LIBS_DEPENDENCIES openblas)
	add_definitions(-DUSE_OPENBLAS)
endif()
if(LINK_LAPACK)
	list(APPEND LIBS_DEPENDENCIES lapack)
endif()

# OpenMP is optional (does not affect OpenMP usage in Armadillo)
if(ENABLE_OPENMP)
	# The -fopenmp flag has be to be inserted after linking Armadillo, so add it to the end of the dependencies list
	list(APPEND LIBS_DEPENDENCIES -fopenmp)
	add_definitions(-DUSE_OPENMP)
endif()

# -------------------------------------------------
# Build targets
# -------------------------------------------------
# Core functionality libraries
add_library(molspin_core STATIC ${SOURCES_TARGET_LIBRARY_MOLSPIN_CORE})
# MolSpin application
add_executable(molspin ${SOURCES_TARGET_APPLICATION_MOLSPIN})
target_link_libraries(molspin ${LIBS_MOLSPIN_CORE} ${LIBS_DEPENDENCIES})
# Tests application
add_executable(tests ${SOURCES_TARGET_APPLICATION_TESTS})
target_link_libraries(tests ${LIBS_MOLSPIN_CORE} ${LIBS_DEPENDENCIES})
