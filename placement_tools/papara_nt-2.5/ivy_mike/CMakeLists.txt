cmake_minimum_required( VERSION 2.8 )

IF(WIN32)
IF(NOT BOOST_ROOT)
set(BOOST_ROOT C:\\src\\boost_1_47_0)
ENDIF()


# include_directories( ${BOOST_ROOT} )
file(GLOB IM_HEADERS src/ivymike/*.h)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_SCL_SECURE_NO_WARNINGS") # as long as there is no support for std::array, these warnings are plain stupid!
ELSE()

  if( USE_CPP11 )
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -march=native")
  else()
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++98 -pedantic -Wall -march=native")
  endif()

ENDIF()

IF( BOOST_ROOT )
include_directories( ${BOOST_ROOT} )
ENDIF()

PROJECT(ivymike)
SET( IVYMIKE_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/src CACHE STRING "IvyMike include dir" FORCE)

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

#ADD_LIBRARY(ivymike src/main.cpp src/LargePhylip.cpp src/time.cpp )
ADD_LIBRARY(ivymike STATIC src/time.cpp src/tree_parser.cpp src/multiple_alignment.cpp src/getopt.cpp src/demangle.cpp src/sdf.cpp src/tree_split_utils.cpp src/LargePhylip.cpp src/large_phylip.cpp ${IM_HEADERS})
set_property(TARGET ivymike PROPERTY CXX_STANDARD 11)

#INSTALL(TARGETS ivymike
#  LIBRARY DESTINATION lib
#  ARCHIVE DESTINATION lib
#)

# ADD_EXECUTABLE( treeparser_test  src/tree_parser.cpp src/tree_parser_test.cpp )
# ADD_EXECUTABLE( sdf_test src/sdf.cpp )


#INSTALL(FILES ivymike/LargePhylip.h
#  DESTINATION include)

