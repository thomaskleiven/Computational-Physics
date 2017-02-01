cmake_minimum_required (VERSION 2.6)
project (oving_1)
find_package(Armadillo REQUIRED)
link_directories("/usr/lib64")

include_directories(${ARMADILLO_INCLUDE_DIRS})
set(lib ${ARMADILLO_LIBRARIES})

find_library( LAPACK lapack )
if ( NOT LAPACK )
  MESSAGE( FATAL_ERROR "Did not find lapack!")
endif()
MESSAGE( STATUS "Found lapack ${LAPACK}")

set( lib ${lib} ${LAPACK} )

add_executable(main main.cpp)
target_link_libraries(main ${lib})
