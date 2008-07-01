cmake_minimum_required(VERSION 2.8)
project(DTRA Fortran C CXX)

#include(FortranCInterface)
#FortranCInterface_HEADER(DTRA/fc.h SYMBOLS mars3d_RPICFD)
#FortranCInterface_VERIFY(CXX)

#include_directories(${DTRA_BINARY_DIR})
#add_executable(mars mars.cpp main.f90 DTRA/fc.h)
#set_property(TARGET mars PROPERTY LINKER_LANGUAGE Fortran)

add_executable(mars2 mars.c main.f90 mars.cpp)
#set_property(TARGET mars2 PROPERTY LINKER_LANGUAGE CXX)
set_property(TARGET mars2 PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(mars2
  #${DTRA_SOURCE_DIR}/../test-sphere/libmars.a
  #${DTRA_SOURCE_DIR}/lib/libmars.a
  -lguide
  -lpthread
  -lz
 # -cxxlib
  )
