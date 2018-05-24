#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(GenFind CXX)

# Add library source directory name that appear in ./src
set(GenFind_LIB_NAMES
  GenFind
  )

# Add executable source directory name that appear in ./src
set(GenFind_EXE_NAMES
  #GenFind
  runtest
  )

# The version number.
set(GenFind_MAJOR_VERSION 0)
set(GenFind_MINOR_VERSION 1)
set(GenFind_PATCH_VERSION 0)
SET(GenFind_VERSION "${GenFind_MAJOR_VERSION}.${GenFind_MINOR_VERSION}.${GenFind_PATCH_VERSION}")

# ---------------------------------------------------------------------------
# load additional OS dependent settings
include(cmake/os.cmake)

# ---------------------------------------------------------------------------
# Configure the find_package outputs for this library
# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR     lib     CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR     bin     CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
   set(DEF_INSTALL_CMAKE_DIR CMake)
else()
   set(DEF_INSTALL_CMAKE_DIR lib/${PROJECT_NAME})
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake
   ${PROJECT_SOURCE_DIR}/cmake/Modules 
   )

if(NOT CMAKE_INSTALL_PREFIX)
   SET(GenFind_PREFIX "/usr/local")
else() 
   SET(GenFind_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

# ---------------------------------------------------------------------------
#
set(GenFind_CXX_FLAGS  " " CACHE PATH "CXX compile flags")
set(GenFind_CXX_LDFLAGS "-L${CMAKE_INSTALL_PREFIX}/lib ${CMAKE_SHARED_LINKER_FLAGS}" CACHE PATH "Linking flags")
set(GenFind_CXX_INC_DIR "-I${CMAKE_INSTALL_PREFIX}/include/GenFind "          CACHE PATH "include dir")
set(GenFind_CXX_LIBS    "" CACHE PATH "libraries")
set(GenFind_FC_FLAGS   "" CACHE PATH "fortran compile flags")
set(GenFind_DATA_DIR $ENV{JLAB_ROOT}/noarch/data CACHE PATH "data directory")

#set (GenFind_COMPILER_FLAGS "-Wall")
set(GenFind_COMPILER_FLAGS   "-Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${GenFind_COMPILER_FLAGS} ")
set(GenFind_PDF_GRID_DIR     "${PROJECT_SOURCE_DIR}/grid")
set(GenFind_GRID_DATA_DIR    "${GenFind_GRID_DATA_DIR} ${GenFind_PDF_GRID_DIR}")
set(CMAKE_USE_ROOT6        ON)

# ---------------------------------------------------------------------------
# Must be out of source build
macro(GenFind_CHECK_OUT_OF_SOURCE_BUILD)
  get_filename_component(bindir_parent ${CMAKE_BINARY_DIR} PATH)
  if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    file(REMOVE_RECURSE ${CMAKE_SOURCE_DIR}/Testing)
    file(REMOVE ${CMAKE_SOURCE_DIR}/DartConfiguration.tcl)
    message(FATAL_ERROR "GenFind should be built as an out of source build, to keep the source directory clean. Please create a extra build directory and run the command 'cmake <path_to_source_dir>' in this newly created directory.  You have also to delete the directory CMakeFiles and the file CMakeCache.txt in the source directory. Otherwise cmake will complain even if you run it from an out-of-source directory.")
  elseif(IS_SYMLINK ${CMAKE_BINARY_DIR} AND CMAKE_SOURCE_DIR STREQUAL bindir_parent)
    message(FATAL_ERROR "GenFind cannot be built from a sub-directory of the source tree that is a symlink. This is a current limitation of CMake. Please create a real build directory and run the command 'cmake <path_to_source_dir>' in this newly created directory.")
  endif()
endmacro()
# NOTE: Doesn't work with older cmake versions:
GenFind_CHECK_OUT_OF_SOURCE_BUILD()

# Check the compiler and set the compile and link flags
#set(CMAKE_BUILD_TYPE Debug)

set(LIBRARY_OUTPUT_PATH    "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
include_directories(APPEND ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Find Libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist MathMore Graf Graf3d
   Gpad GX11 GeomPainter X3d Gviz3d RGL
   Eve EG Foam GenVector Geom GLEW Gui HistPainter MathCore Matrix MLP Net
   Physics Spectrum Thread Tree TreeViewer)
include(${ROOT_USE_FILE})
include(cmake/root.cmake)

# ----------
# CLHEP
#set (CLHEP_CONFIG_DEBUG 1)
#find_package(CLHEP  REQUIRED)
##  CLHEP_INCLUDE_DIRS - include directories for CLHEP
##  CLHEP_DEFINITIONS  - compile definitions needed to use CLHEP
##  CLHEP_LIBRARIES    - libraries to link against
#set(CLHEP_LIBRARY_DIR $ENV{JLAB_PREFIX}/lib)
#include_directories(${CLHEP_INCLUDE_DIR})
##message(" CLHEP_DEFINITIONS = ${CLHEP_DEFINITIONS} " )

# ----------
## EVIO
##set(EVIO_INCLUDE_DIR $ENV{JLAB_PREFIX}/include)
##set(EVIO_LIBRARY_DIR $ENV{JLAB_PREFIX}/lib)
##set(EVIO_LIBRARIES  evio evioxx expat)

set(LINK_LIBRARIES 
   ${LINK_LIBRARIES}
   ${ROOT_LIBRARIES}
   #${EVIO_LIBRARIES}
   #${CLHEP_LIBRARIES}
   )

#----------------------------------------------------------------------------
# not sure this is the best way to do this....
set(GenFind_FC_FLAGS "${GenFind_FC_FLAGS} ${CMAKE_Fortran_FLAGS}")
#set(CMAKE_Fortran_FLAGS " -D\"GenFind_PDF_GRID_DIR='${GenFind_PDF_GRID_DIR}'\" -g -O2 -Wall -fPIC -Iinclude -ffixed-line-length-none -x f77-cpp-input ")
set(CMAKE_Fortran_FLAGS " -D\"GenFind_PDF_GRID_DIR='${GenFind_PDF_GRID_DIR}'\" -g -O2 -Wall  -ffixed-line-length-none -x f77-cpp-input ")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}")
set(CMAKE_Fortran_FLAGS " -D\"GenFind_PDF_GRID_DIR='${GenFind_PDF_GRID_DIR}'\" -g -O2 -Wall  -ffixed-line-length-none -x f77-cpp-input ")
set(GenFind_CXX_LDFLAGS   "${GenFind_CXX_LDFLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")

# -------------------------------------------------------
# Make string of libs for commandline config tool
set(GenFind_SHARED_LIB_NAMES)
foreach(aSharedLib ${GenFind_LIB_NAMES})
   set(GenFind_CXX_LIBS          "${GenFind_CXX_LIBS} -l${aSharedLib} ")
   set(GenFind_SHARED_LIB_NAMES   ${GenFind_SHARED_LIB_NAMES} ${aSharedLib})
endforeach(aSharedLib)

set(GenFind_SHARED_LIBS)
set(GenFind_LIBRARIES)
set(GenFind_ROOT_DICTS)
set(GenFind_PCM_FILES)

#---------------------------------------
foreach(aSharedLib ${GenFind_LIB_NAMES})
   set(GenFind_LIBRARIES ${GenFind_LIBRARIES} ${aSharedLib})
   #set(GenFind_SHARED_LIBS ${GenFind_SHARED_LIBS} "${CMAKE_BINARY_DIR}/lib/GenFind${aSharedLib}.so")
   set(GenFind_ROOT_DICTS ${GenFind_ROOT_DICTS} "GenFind${aSharedLib}Dict.h")
   set(GenFind_ROOT_DICT_SRCS ${GenFind_ROOT_DICT_SRCS} "GenFind${aSharedLib}Dict.cxx")
   set(GenFind_PCM_FILES ${GenFind_PCM_FILES} "GenFind${aSharedLib}Dict_rdict.pcm")
   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/${aSharedLib}/include)
endforeach(aSharedLib)

# -------------------------------------------------------
# configure a header file to pass some of the CMake settings to the source code
configure_file (
   "${PROJECT_SOURCE_DIR}/src/GenFindConfig.h.in"
   "${PROJECT_BINARY_DIR}/GenFindConfig.h"
   )

#----------------------------------------------------------------------------
# Add the config tool directory
#add_executable(GenFind-config GenFind-config.cxx )
add_subdirectory(src/config)

#-----------------------------------------
# add the library sub directories
foreach(aSharedLib ${GenFind_LIB_NAMES})
   #include(${aSharedLib})
   add_subdirectory("src/${aSharedLib}")
endforeach(aSharedLib)

# add the executable sub directories
foreach(anExeName ${GenFind_EXE_NAMES})
  #include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/${anExeName}/include)
  add_subdirectory("src/${anExeName}")
endforeach(anExeName)

FILE(GLOB GENERATED_HEADERS 
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h" 
  "${PROJECT_BINARY_DIR}/GenFindConfig.h")
FILE(GLOB PCM_FILES "${PROJECT_BINARY_DIR}/*.pcm")

# -----------------------------------------------------------
# Set the library version in the main CMakeLists.txt
SET(GenFind_LIBRARY_PROPERTIES ${GenFind_LIBRARY_PROPERTIES}
  VERSION "${GenFind_VERSION}"
  SOVERSION "${GenFind_MAJOR_VERSION}"
  SUFFIX ".so"
  )

#----------------------------------------------------------------------------
# Configure and install "share" files.  
# These include mostly root scripts and examples.
set(SHARE_SCRIPTS
  #examples/example1.cxx
  #share/plot_field.cxx
  #share/plot_field_xy.cxx
  )
foreach(_script ${SHARE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
install(FILES ${SHARE_SCRIPTS} DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}" )

# install examples
set(GenFind_EXAMPLES 
  examples/example1.cxx 
  )
install(FILES ${GenFind_EXAMPLES} DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/examples" )

# ----------------------------------------------------------------------------
# The library's cmake files
# Add all targets to the build-tree export set
#export(TARGETS GenFind FILE "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")
# TODO: replace "GenFind" with smaller libraries

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE ${PROJECT_NAME})

set(CONF_LIBRARIES 
   ${GenFind_SHARED_LIB_NAMES}  # TODO: Break big library into smaller ones
   )

set(CONF_EXECUTABLES 
   GenFind-config
   )

set(CONF_INCLUDE_DIRS 
   ${CMAKE_INSTALL_PREFIX}/include
   ${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}
   )

set(CONF_LIBRARY_DIRS 
   ${CMAKE_INSTALL_PREFIX}/lib
   )

configure_file(src/config/${PROJECT_NAME}Config.cmake.in 
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake" @ONLY)
configure_file(src/config/${PROJECT_NAME}ConfigVersion.cmake.in 
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}ConfigVersion.cmake" @ONLY)
configure_file(src/config/${PROJECT_NAME}UseFile.cmake.in 
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}UseFile.cmake" @ONLY)

# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}ConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}UseFile.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}Targets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)

#------------------------------------
# Debug
MESSAGE("          GenFind_CXX_LIBS: ${GenFind_CXX_LIBS}")
MESSAGE("       GenFind_CXX_LDFLAGS: ${GenFind_CXX_LDFLAGS}")
MESSAGE("       GenFind_CXX_INC_DIR: ${GenFind_CXX_INC_DIR}")
MESSAGE("             CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
MESSAGE("              FC_COMPILER: ${CMAKE_FC_COMPILER}")
MESSAGE("      CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
MESSAGE("                CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
MESSAGE("         CMAKE_LINK_FLAGS: ${CMAKE_LINK_FLAGS}")
MESSAGE("CMAKE_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}")

message("GenFind_MAJOR_VERSION   ${GenFind_MAJOR_VERSION}")
message("GenFind_MINOR_VERSION   ${GenFind_MINOR_VERSION}")
message("GenFind_PATCH_VERSION   ${GenFind_PATCH_VERSION}")
message("GenFind_DATA_DIR      ${GenFind_DATA_DIR}")
message("GenFind_CXX_FLAGS    ${GenFind_CXX_FLAGS}")
message("GenFind_CXX_LDFLAGS   ${GenFind_CXX_LDFLAGS}")
message("GenFind_CXX_LIBS      ${GenFind_CXX_LIBS}")
message("GenFind_CXX_INC_DIR   ${GenFind_CXX_INC_DIR}")
message("GenFind_PREFIX        ${CMAKE_INSTALL_PREFIX}")

