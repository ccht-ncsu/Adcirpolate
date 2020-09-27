# What is adcirpolate ?
*adcirpolate* is a module for interpolating hotstart files between different ADCIRC meshes.
Hence, we can switch between different meshes to avoid unnecessary computation on meshes with
high resolution in areas where no important phenomenon is going on. If you do not know about
ADCIRC, or its hotstart files, you are in the wrong place !

- [Prerequisites](##Prerequisites)
- [CMake and Fortran version requirements](##CMake-and-Fortran-version-requirements)
- [Running the code](##Running-the-code)

## Prerequisites
In order to run this code, you need MPI and ESMF libraries. The program is to be built with the same MPI library that is used to build ESMF. You also need to set the following environment variables, and run cmake:

    export FC=mpif90
    export ESMF_CONFIG_FILE=/path/to/esmf.mk
    cmake /path/to/CMakeLists.txt
    make all

## CMake and Fortran version requirements
The current version of Adcirpolate has been tested on multiple HPC clusters with different versions of CMake and GNU and Intel Fortran compilers. Understandably, some tunings to the CMakeLists.txt might be required for clusters with less up-to-date software. In case you are having trouble compiling the code on your machine, please submit an issue with output log from CMake or Fortran compiler.

## Running the code
Adcirpolate should be executed in a directory with two folders named `coarse` and `fine`. Typically, we transfer the results from hotstart file of the coarse mesh to a hotstart file of the fine mesh. The two folders contain the partitioned ADCIRC meshes. In the `coarse` directory, the hotstart file of the coarse mesh exists. After running, the program creates a hotstart file for the fine mesh.
