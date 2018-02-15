cmake_minimum_required (VERSION 2.8)
project (adcirpolate Fortran)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (ESMF)

if (DEFINED WITH_HDF5 AND NOT DEFINED ENV{HDF5_DIR})
  message (FATAL_ERROR "\n${Esc}[1;31m!! Error: You have enabled HDF5 linking, \
                   but you have not defined HDF5_DIR. Try something like: ${Esc}[m\
                   export HDF_DIR=/path/to/hdf5 && cmake ${CMAKE_SOURCE_DIR}")
endif ()

if (DEFINED ENV{HDF5_DIR})
  set (HDF5_ROOT ${HDF5_DIR})
  find_package (HDF5)
endif ()

add_executable(${PROJECT_NAME} adcirpolate.f90
               main.f90 wetdry.f90 errors_and_msgs.f90
               c_preprocessor.h)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_COMPILER_LINE}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_LINK_LINE} -g -cpp")

#
# Adding record length interpretation for Intel
#
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                            -assume byterecl")
  message ("+> Compiler is intel and -assume byterecl is used.")
endif ()

#
# Building in debug mode
#
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG_MODE)
  message("Building the software in debug mode.")
endif()

add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT " ")

