cmake_minimum_required (VERSION 2.8.7 FATAL_ERROR)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)

# main compile alternators and their defaults
set (ENABLE_MPI       OFF)
set (ENABLE_OMP       ON)
set (ENABLE_PREFETCH  ON)

if ( $ENV{EPA_HYBRID} )
	set (ENABLE_MPI       ON)
	set (ENABLE_OMP       ON)
	set (ENABLE_PREFETCH  ON)
endif ()

if ( $ENV{EPA_SERIAL} )
    set (ENABLE_MPI       OFF)
    set (ENABLE_OMP       OFF)
    set (ENABLE_PREFETCH  OFF)
endif ()

project ( epa CXX C )

set (epa_VERSION_MAJOR 0)
set (epa_VERSION_MINOR 2)
set (epa_VERSION_PATCH 1)
set (epa_VERSION_BETA "beta")

set (epa_PROJECT_VERSION "${epa_VERSION_MAJOR}")
set (epa_PROJECT_VERSION "${epa_PROJECT_VERSION}.${epa_VERSION_MINOR}")
set (epa_PROJECT_VERSION "${epa_PROJECT_VERSION}.${epa_VERSION_PATCH}")
set (epa_PROJECT_VERSION "${epa_PROJECT_VERSION}-${epa_VERSION_BETA}")

#check for minimum compiler versions
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(EPA_COMPILER_TARGET_VERSION "4.9.2")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(EPA_COMPILER_TARGET_VERSION "3.8")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS EPA_COMPILER_TARGET_VERSION)
    message (FATAL_ERROR "${CMAKE_CXX_COMPILER_ID} compiler too old! Minimum required: ${EPA_COMPILER_TARGET_VERSION}")
  endif()
endif()

message(STATUS "Building EPA-ng version: ${epa_PROJECT_VERSION}")

set (CMAKE_BUILD_TYPE RELEASE)

if ( $ENV{EPA_DEBUG} )
  set (CMAKE_BUILD_TYPE DEBUG)
endif ()

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set (WARN_FLAGS               "-Wall -Wextra")
set (CMAKE_CXX_FLAGS          "-std=c++14 ${WARN_FLAGS} -DEPA_VERSION='\"${epa_PROJECT_VERSION}\"'")

set (CMAKE_CXX_FLAGS_DEBUG    "-O0 -g -ggdb3 -DDEBUG -D_GLIBCXX_DEBUG")
set (CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

if(ENABLE_PREFETCH)
  set (CMAKE_THREAD_PREFER_PTHREAD ON)
  set (THREADS_PREFER_PTHREAD_FLAG ON)
  find_package (Threads)
  message(STATUS "Enabling Prefetching")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__PREFETCH")
endif()

if(ENABLE_OMP)
  message(STATUS "Checking for OpenMP")
  find_package( OpenMP QUIET )
  if(OPENMP_FOUND)
    message(STATUS "Checking for OpenMP -- found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D__OMP")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else()
    message(STATUS "Checking for OpenMP -- NOT FOUND")
  endif()
endif()

if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_CXX_FOUND)
    set( ENV{OMPI_CXX} ${CMAKE_CXX_COMPILER} PARENT_SCOPE )
    #set( CMAKE_CXX_COMPILER       "mpicxx")
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILER_FLAGS} -D__MPI")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
  endif()
endif()

# =============================================
#   Download Dependencies
# =============================================

# Downloads dependencies if they are not there. We use a special script for this.
include( "${CMAKE_CURRENT_LIST_DIR}/libs/cmake/DownloadDependency.cmake" )

# Store the commit hashes of our dependencies.
# These are replaced by libs/cmake/update_dependencies.sh to the hashes that are currently checked out.
# Thus, do not replace the hashes manually!
SET( CLI11_COMMIT_HASH "d465414455e01c9e0753efbbf5bde16685bc7b90" ) #CLI11_COMMIT_HASH#
SET( genesis_COMMIT_HASH "fcfa3e2d45df6e4b64427c26cdf3fc4e2767c429" ) #genesis_COMMIT_HASH#
SET( pll-modules_COMMIT_HASH "f6174cdff472be5a43220e1581d5d4cb7d4ef18c" ) #pll-modules_COMMIT_HASH#
SET( libpll_COMMIT_HASH "f89f1f40842ecba4ed0b1e6a87a688f9a46865fb" ) #libpll_COMMIT_HASH#
SET( cereal_COMMIT_HASH "51cbda5f30e56c801c07fe3d3aba5d7fb9e6cca4" ) #cereal_COMMIT_HASH#

# Call the github download function, which takes four arguments:
#  - LIBPATH    : Path to the libracy dir where dependencies are stored.
#  - LIBNAME    : Name of the dependency, that is, the name of its main directory within the ${LIBPATH}.
#  - TESTFILE   : A testfile to check if the dependency is already there.
#  - REPOSITORY : Github repository (user/reponame)
#  - COMMITHASH : Hash of the commit to check out
DOWNLOAD_GITHUB_DEPENDENCY( "${PROJECT_SOURCE_DIR}/libs" "CLI11"      "CMakeLists.txt"  "CLIUtils/CLI11"        ${CLI11_COMMIT_HASH} )
DOWNLOAD_GITHUB_DEPENDENCY( "${PROJECT_SOURCE_DIR}/libs" "genesis"    "CMakeLists.txt"  "lczech/genesis"        ${genesis_COMMIT_HASH} )
DOWNLOAD_GITHUB_DEPENDENCY( "${PROJECT_SOURCE_DIR}/libs" "pll-modules""README.md"       "ddarriba/pll-modules"  ${pll-modules_COMMIT_HASH} )
DOWNLOAD_GITHUB_DEPENDENCY( "${PROJECT_SOURCE_DIR}/libs/pll-modules/libs"
                                                         "libpll"     "README.md"       "xflouris/libpll-2"     ${libpll_COMMIT_HASH} )
DOWNLOAD_GITHUB_DEPENDENCY( "${PROJECT_SOURCE_DIR}/libs" "cereal"     "README.md"       "USCiLab/cereal"        ${cereal_COMMIT_HASH} )

# =============================================
#         Build pll-modules
# =============================================

add_subdirectory(${PROJECT_SOURCE_DIR}/libs/pll-modules)

include_directories(${PLLMOD_EXPORTED_INCLUDE_PATH})

# =============================================
#         Build genesis
# =============================================

# pass special info to genesis
if(ENABLE_PREFETCH)
  set(GENESIS_USE_THREADS ON CACHE BOOL "Passing GENESIS_USE_THREADS=ON to genesis")
else()
  set(GENESIS_USE_THREADS OFF CACHE BOOL "Passing GENESIS_USE_THREADS=OFF to genesis")
endif()

set(GENESIS_USE_OPENMP OFF CACHE BOOL "Passing GENESIS_USE_OPENMP=OFF to genesis")

# Add Genesis as dependency
add_subdirectory(${PROJECT_SOURCE_DIR}/libs/genesis)

# Use everything that Genesis exports, just to be sure that we use the same setup.
add_definitions( ${GENESIS_DEFINITIONS} )
include_directories( ${GENESIS_INCLUDE_DIR} )
set( CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${GENESIS_CXX_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GENESIS_EXE_LINKER_FLAGS}" )

message(STATUS "Using flags: ${CMAKE_CXX_FLAGS}")

# =============================================
# Add include directories for other dependencies
# =============================================

include_directories(${PROJECT_SOURCE_DIR}/libs/cereal/include)
include_directories(${PROJECT_SOURCE_DIR}/libs/CLI11/include)

# =============================================
#         Recurse and compile main code
# =============================================

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/test/src)
