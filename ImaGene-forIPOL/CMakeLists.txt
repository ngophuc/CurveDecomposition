# The name of our project is "ImaGene". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${IMAGENE_SOURCE_DIR} and to the root binary directory of the
# project as ${IMAGENE_BINARY_DIR}.

# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
cmake_minimum_required (VERSION 2.4) 
project (ImaGene)

message( "-------------------------------------------------------------------------------" )
message( "ImaGene: written by Jacques-Olivier Lachaud (Univ. of Savoie)." )
message( "       You may specify the following options:" )
message( "         - IMAGENE64: 64 bits version of ImaGene." )
message( "         - TRACE_LEVEL = 0 (no) - 3 (verbose): verbosity level for trace.")
message( "         - TIMING_LEVEL = 0 (no) - 3 (all):  timings.")
message( "         - CMAKE_BUILD_TYPE = Debug|Release: choose configuration.")
message( "         - CMAKE_INSTALL_PREFIX = <path>: the place to install ImaGene (/usr/local).")
message( "       with the syntax (for instance)")  
message( "         cmake -DCMAKE_BUILD_TYPE:STRING=Debug")
message( "-------------------------------------------------------------------------------" )
message("Source DIR is ${PROJECT_SOURCE_DIR}")
message("Binary DIR is ${PROJECT_BINARY_DIR}")

#add_definitions("-DTRACE_DMLP")

set(BUILD_SHARED_LIBS true)
if (DEFINED TRACE_LEVEL)
else (DEFINED TRACE_LEVEL)
   set(TRACE_LEVEL 0)
endif (DEFINED TRACE_LEVEL)
if (DEFINED TIMING_LEVEL)
else (DEFINED TIMING_LEVEL)
     set(TIMING_LEVEL 0)
endif (DEFINED TIMING_LEVEL)

if (DEFINED IMAGENE64)
   add_definitions( -DIMAGENE64 )
   message( STATUS "64 bits mode chosen." )
   set( SUFFIXLIB "64" )
   set( SUFFIXBIN "64" )
else (DEFINED IMAGENE64)
   message( STATUS "32 bits mode chosen." )
   set( SUFFIXLIB "" )
   set( SUFFIXBIN "" )
endif (DEFINED IMAGENE64)


set(INSTALL_DOC_PATH ${CMAKE_INSTALL_PREFIX}/doc/${CMAKE_PROJECT_NAME} )
if (${CMAKE_BUILD_TYPE} MATCHES D[eE][bB][uU][gG])
   set(CMAKE_BUILD_TYPE Debug)
   set( CMAKE_C_FLAGS_DEBUG "-g -DDEBUG" )
   set( CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG" )
   set( SUFFIXLIB "${SUFFIXLIB}-debug" )
   set( SUFFIXBIN "${SUFFIXBIN}-debug" )
   message( STATUS "DEBUG mode with flags = ${CMAKE_CXX_FLAGS_DEBUG}" )
else (${CMAKE_BUILD_TYPE} MATCHES D[eE][bB][uU][gG])
   set(CMAKE_BUILD_TYPE Release)
   set( CMAKE_C_FLAGS_RELEASE "-O3 -UDEBUG -DNO_DEBUG" )
   set( CMAKE_CXX_FLAGS_RELEASE "-O3 -UDEBUG -DNO_DEBUG" )
   message( STATUS "RELEASE mode with flags = ${CMAKE_CXX_FLAGS_RELEASE}" )
endif (${CMAKE_BUILD_TYPE} MATCHES D[eE][bB][uU][gG])
message( STATUS "Libs extension is '${SUFFIXLIB}'. Bins extension is '${SUFFIXBIN}'." )
# -----------------------------------------------------------------------------
# Defines some variables
# -----------------------------------------------------------------------------
set ( LIBIMAGENE_NAME "ImaGene${SUFFIXLIB}" )

# add_definitions( "-DBIG_INTEGER_IS_LONGLONG" )

# -----------------------------------------------------------------------------
# Display some information
# -----------------------------------------------------------------------------
message(STATUS "Installation prefix directory is " ${CMAKE_INSTALL_PREFIX})
message(STATUS "Host system is " ${CMAKE_HOST_SYSTEM} " with processor " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "Target system is " ${CMAKE_SYSTEM} " with processor " ${CMAKE_SYSTEM_PROCESSOR})
# message(STATUS "Build type is " ${CMAKE_BUILD_TYPE})

# -----------------------------------------------------------------------------
# Check packages
# -----------------------------------------------------------------------------
# Look for modules
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
list( APPEND CMAKE_MODULE_PATH . )
#message("CMAKE_MODULE_PATH=" ${CMAKE_MODULE_PATH})





# nead file 'doxy.config.in' for documentation
include(doc/doxygen.cmake)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
INCLUDE("${CMAKE_MODULE_PATH}/TargetDoxygenDoc.cmake" OPTIONAL)
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/TargetUninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/TargetUninstall.cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/TargetUninstall.cmake")

# find_package(DOXYGEN)
# set(DOXYGEN_CONFIG_FILE doc/Doxyfile)
# if (DOXYGEN_FOUND)
#    SET(DOXYFILE_FOUND false)
#    IF (EXISTS ${PROJECT_SOURCE_DIR}/${DOXYGEN_CONFIG_FILE})
#       SET(DOXYFILE_FOUND true)
#    ENDIF(EXISTS ${PROJECT_SOURCE_DIR}/${DOXYGEN_CONFIG_FILE})
#    IF (DOXYFILE_FOUND)
#       ADD_CUSTOM_TARGET(doc ${DOXYGEN_EXECUTABLE} "${PROJECT_SOURCE_DIR}/${DOXYGEN_CONFIG_FILE}" COMMENT "Create Doxygen documentation")
#       set_target_properties(doc PROPERTIES EXCLUDE_FROM_ALL 1)
#    endif (DOXYFILE_FOUND)
# endif (DOXYGEN_FOUND)

#------------------------------------------------------------------------------
# Add subdirectories
#------------------------------------------------------------------------------
add_subdirectory (src)
add_subdirectory (include)
add_subdirectory (bin)
add_subdirectory (examples)
add_subdirectory (demoIPOL)

# add_subdirectory (digitalPlaneRecognition)
# add_subdirectory (doc)

#------------------------------------------------------------------------------
# CTest Stuff
#------------------------------------------------------------------------------
set(VALGRIND_COMMAND_OPTIONS "-q --tool=memcheck --leak-check=yes --show-reachable=yes --workaround-gcc296-bugs=yes --num-callers=50")
ENABLE_TESTING()
add_subdirectory (tests)

#Following commands are used to enable the gcc coverage analysis in CTEST (Debug build type only)
if (CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-parameter -Wunused-function        -Wunused -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
  SET(CMAKE_LD_FLAGS_DEBUG "${CMAKE_LD_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif (CMAKE_COMPILER_IS_GNUCXX)
if (CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
  SET(CMAKE_LD_FLAGS_DEBUG "${CMAKE_LD_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
endif (CMAKE_COMPILER_IS_GNUCC)


#------------------------------------------------------------------------------
# Prepare installations concerning global scope.
#------------------------------------------------------------------------------
if (CMAKE_MAJOR_VERSION EQUAL 2)
   if (CMAKE_MINOR_VERSION GREATER 5)
      install(DIRECTORY include/ DESTINATION include 
      			FILES_MATCHING PATTERN "*.*h")
   else (CMAKE_MINOR_VERSION GREATER 5)
      install(DIRECTORY include/ DESTINATION include 
      			PATTERN "*.*h")
   endif (CMAKE_MINOR_VERSION GREATER 5)
endif (CMAKE_MAJOR_VERSION EQUAL 2)


#------------------------------------------------------------------------------
# Some directories and files should also be cleaned when invoking 'make clean'
#------------------------------------------------------------------------------
ADD_CUSTOM_TARGET(distclean
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/TargetDistclean.cmake")


#------------------------------------------------------------------------------
# Packaging
#------------------------------------------------------------------------------
INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ImaGene: a library for manipulating digital spaces.")
SET(CPACK_PACKAGE_VENDOR "Jacques-Olivier Lachaud")
# SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
# SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYRIGHT")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
# SET(CPACK_PACKAGE_EXECUTABLES "MyExecutable" "My Executable")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES _CPack_Packages;.*.tar.gz;.*.tgz;.*.tar.Z;.*.dmg;[.]svn;CMakeCache.txt;${PROJECT_NAME}.build;Debug;Release)

INCLUDE(CPack)
