cmake_minimum_required (VERSION 2.8)
project (adcirpolate Fortran)

# Defining the ${Esc} for syntax coloring.
string(ASCII 27 Esc)

# Checking if ESMF exists
if (NOT DEFINED ENV{ESMF_CONFIG_FILE} AND NOT DEFINED ESMF_CONFIG_FILE)
  message (FATAL_ERROR "\n${Esc}[1;31m!! Error: You need ESMF library to \
                   run this program. please set the environment \
                   variable ESMF_CONFIG_FILE to point to esmf.mk in \
                   your ESMF installation directory. \
                   Try something like: ${Esc}[m\
                   export ESMF_CONFIG_FILE=/path/to/esmf.config && cmake ${CMAKE_SOURCE_DIR}")
endif ()

if (NOT EXISTS $ENV{ESMF_CONFIG_FILE} AND NOT EXISTS ${ESMF_CONFIG_FILE})
  message (FATAL_ERROR "${Esc}[1;31m Error: esmf.mk file is not found at \
                        ${ESMF_CONFIG_FILE} ${Esc}[m")
else ()
  message ("+>${Esc}[1;32m The config file for ESMF library is found.${Esc}[m")
endif ()

if (DEFINED ENV{ESMF_CONFIG_FILE})
  set(ESMF_CONFIG_FILE $ENV{ESMF_CONFIG_FILE} CACHE STRING "")
endif ()
set(ESMF_CONFIG_FILE ${ESMF_CONFIG_FILE} CACHE STRING "")

file(STRINGS "${ESMF_CONFIG_FILE}" all_vars)
foreach(str ${all_vars})
  string(REGEX MATCH "^[^#]" def ${str})
  if (def)
    string(REGEX MATCH "^[^=]+" var_name ${str})
    string(REGEX MATCH "=(.+)$" var_def ${str})
    set(var_def ${CMAKE_MATCH_1})
    set(${var_name} ${var_def})
  endif()
endforeach()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions(-DDEBUG_MODE)
  message("Building the software in debug mode.")
endif()

add_executable(${PROJECT_NAME} main.f90)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_F90COMPILEOPTS}\
                         ${ESMF_F90COMPILEPATHS}\
                         ${ESMF_F90COMPILEFREENOCPP}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_F90LINKOPTS}\
                         ${ESMF_F90LINKRPATHS}\
                         ${ESMF_F90LINKPATHS}\
                         ${ESMF_F90ESMFLINKLIBS} -g -cpp")

add_custom_target(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT " ")

