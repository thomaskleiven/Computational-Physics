
cmake_minimum_required (VERSION 2.6)
project (Task_1)
find_package(Armadillo REQUIRED)
link_directories("/usr/lib64")

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set (CMAKE -DCMAKE_BUILD_TYPE=RELEASE)

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

add_executable(main Main.cpp Activator.cpp Lattice.cpp)
target_link_libraries(main ${LIB})
