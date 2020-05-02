cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(mapleLRB)

#--------------------------------------------------------------------------------------------------
# Configurable options:

option(STATIC_BINARIES "Link binaries statically." ON)
option(USE_SORELEASE   "Use SORELEASE in shared library filename." ON)

#--------------------------------------------------------------------------------------------------
# Library version:

set(MAPLE_SOMAJOR   2)
set(MAPLE_SOMINOR   1)
set(MAPLE_SORELEASE 0)

# Compute VERSION and SOVERSION:
if (USE_SORELEASE)
  set(MAPLE_VERSION ${MAPLE_SOMAJOR}.${MAPLE_SOMINOR}.${MAPLE_SORELEASE})
else()
  set(MAPLE_VERSION ${MAPLE_SOMAJOR}.${MAPLE_SOMINOR})
endif()
set(MAPLE_SOVERSION ${MAPLE_SOMAJOR})

#--------------------------------------------------------------------------------------------------
# Dependencies:

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIR})
include_directories(${mapleLRB_SOURCE_DIR})

#--------------------------------------------------------------------------------------------------
# Compile flags:

add_definitions(-D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS)

#--------------------------------------------------------------------------------------------------
# Build Targets:

set(MAPLE_LIB_SOURCES
    mapleLRB/utils/Options.cc
    mapleLRB/utils/System.cc
    mapleLRB/core/Solver.cc
    mapleLRB/simp/SimpSolver.cc)

add_library(mapleLRB-lib-static STATIC ${MAPLE_LIB_SOURCES})
add_library(mapleLRB-lib-shared SHARED ${MAPLE_LIB_SOURCES})

target_link_libraries(mapleLRB-lib-shared ${ZLIB_LIBRARY})
target_link_libraries(mapleLRB-lib-static ${ZLIB_LIBRARY})

add_executable(mapleLRB_core mapleLRB/core/Main.cc)
add_executable(mapleLRB_simp mapleLRB/simp/Main.cc)

if(STATIC_BINARIES)
  target_link_libraries(mapleLRB_core mapleLRB-lib-static)
  target_link_libraries(mapleLRB_simp mapleLRB-lib-static)
else()
  target_link_libraries(mapleLRB_core mapleLRB-lib-shared)
  target_link_libraries(mapleLRB_simp mapleLRB-lib-shared)
endif()

set_target_properties(mapleLRB-lib-static PROPERTIES OUTPUT_NAME "mapleLRB")
set_target_properties(mapleLRB-lib-shared
  PROPERTIES
    OUTPUT_NAME "mapleLRB" 
    VERSION ${MAPLE_VERSION}
    SOVERSION ${MAPLE_SOVERSION})

set_target_properties(mapleLRB_simp       PROPERTIES OUTPUT_NAME "mapleLRB")

#--------------------------------------------------------------------------------------------------
# Installation targets:

install(TARGETS mapleLRB-lib-static mapleLRB-lib-shared mapleLRB_core mapleLRB_simp 
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(DIRECTORY mapleLRB/mtl mapleLRB/utils mapleLRB/core mapleLRB/simp
        DESTINATION include/mapleLRB
        FILES_MATCHING PATTERN "*.h")
