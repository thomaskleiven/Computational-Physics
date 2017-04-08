
cmake_minimum_required (VERSION 2.6)
project (src)
find_package(Armadillo REQUIRED)
link_directories("/usr/lib64")

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

include_directories( include )

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

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory( src )
set( LIB ${LIB} oving_3 )

add_executable( main.out bin/Main.cpp )
target_link_libraries( main.out ${LIB} )

add_executable( Step_evolution.out bin/Step_evolution.cpp )
target_link_libraries ( Step_evolution.out ${LIB} )
