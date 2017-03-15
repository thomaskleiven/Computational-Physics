cmake_minimum_required (VERSION 2.6)
project (Task_1)
find_package(Armadillo REQUIRED)
link_directories("/usr/lib64")

include_directories(${ARMADILLO_INCLUDE_DIRS})
set(LIB ${ARMADILLO_LIBRARIES})

find_library( LAPACK lapack )
if ( NOT LAPACK )
  MESSAGE( FATAL_ERROR "Did not find lapack!")
endif()
MESSAGE( STATUS "Found lapack ${LAPACK}")
set( LIB ${LIB} ${LAPACK} )

find_library( GSL gsl )
if( NOT GSL )
  MESSAGE( FATAL_ERROR "Did not find GSL!")
endif()
MESSAGE( STATUS "Found GSL ${GSL}")
set( LIB ${LIB} ${GSL} )

add_executable(main Main.cpp Lattice.cpp)
target_link_libraries(main ${LIB})
