include_directories (${PROJECT_SOURCE_DIR}/src)

file (GLOB_RECURSE epa_sources ${PROJECT_SOURCE_DIR}/src/*.cpp)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable        (epa_module ${epa_sources})

target_link_libraries (epa_module ${GENESIS_LINK_LIBRARIES} )
target_link_libraries (epa_module ${PLLMOD_EXPORTED_LIBS_STATIC})
target_link_libraries (epa_module m)

if(ENABLE_PREFETCH)
  target_link_libraries (epa_module ${CMAKE_THREAD_LIBS_INIT})
endif()

if(ENABLE_MPI)
  if(MPI_CXX_FOUND)
  target_link_libraries (epa_module ${MPI_CXX_LIBRARIES})
  endif()

  if(MPI_COMPILE_FLAGS)
    set_target_properties(epa_module PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  if(MPI_LINK_FLAGS)
    set_target_properties(epa_module PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif()

set_target_properties (epa_module PROPERTIES OUTPUT_NAME epa-ng)
set_target_properties (epa_module PROPERTIES PREFIX "")
