cmake_minimum_required (VERSION 3.0)
project (adcirpolate LANGUAGES Fortran)
enable_language(Fortran)

option(USE_NETCDF "Build with NetCDF" ON)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (ESMF)

add_executable(${PROJECT_NAME} adcirpolate.f90
               main.f90 wetdry.f90 errors_and_msgs.f90
               splitter.f90 netcdfio.f90
               c_preprocessor.h)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_COMPILER_LINE}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}\
                         ${ESMF_LINK_LINE} -g -cpp")

if (${USE_NETCDF})
    set(NETCDF_F90 "YES")
    find_package(NetCDF REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS Fortran)
    target_include_directories(${PROJECT_NAME} PUBLIC ${NETCDF_INCLUDES} ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PUBLIC ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES})
    target_compile_definitions(${PROJECT_NAME} PUBLIC NETCDF)
endif ()

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