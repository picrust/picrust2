find_package (GTest)

if(NOT GTEST_FOUND)
    message (STATUS "GTest not found")
    message (WARNING "Skipping building tests.")
    return()
endif()

message (STATUS "Building tests")

include_directories (${PROJECT_SOURCE_DIR}/src)

file (GLOB_RECURSE epa_test_sources ${PROJECT_SOURCE_DIR}/test/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.cpp)

# sources list now has 2 Main.cpp, old has to be removed
list(REMOVE_ITEM epa_test_sources "${PROJECT_SOURCE_DIR}/src/main.cpp")

include_directories (${PROJECT_SOURCE_DIR})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/test/bin)

add_executable        (epa_test_module ${epa_test_sources})

target_link_libraries (epa_test_module ${GENESIS_LINK_LIBRARIES} )
target_link_libraries (epa_test_module ${PLLMOD_EXPORTED_LIBS_STATIC})
target_link_libraries (epa_test_module m)

target_link_libraries (epa_test_module ${GTEST_BOTH_LIBRARIES})

# if(ENABLE_PREFETCH)
target_link_libraries (epa_test_module ${CMAKE_THREAD_LIBS_INIT})
# endif()

if(ENABLE_MPI)
  if(MPI_CXX_FOUND)
  target_link_libraries (epa_test_module ${MPI_CXX_LIBRARIES})
  endif()

  if(MPI_COMPILE_FLAGS)
    set_target_properties(epa_test_module PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(epa_test_module PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

set_target_properties (epa_test_module PROPERTIES OUTPUT_NAME epa_test)
set_target_properties (epa_test_module PROPERTIES PREFIX "")


add_test (epa_test ${PROJECT_SOURCE_DIR}/test/bin/epa_test)
