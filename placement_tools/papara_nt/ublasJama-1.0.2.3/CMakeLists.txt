IF(WIN32)
IF(NOT BOOST_ROOT)
set(BOOST_ROOT C:\\src\\boost_1_47_0)
ENDIF()

include_directories( ${BOOST_ROOT} )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_SCL_SECURE_NO_WARNINGS") # as long as there is no support for std::array, these warnings are plain stupid!
ELSE()
include_directories( ${BOOST_ROOT} )
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -march=native")

ENDIF()

add_library( ublas_jama CholeskyDecomposition.cpp  EigenvalueDecomposition.cpp  LUDecomposition.cpp  QRDecomposition.cpp )
